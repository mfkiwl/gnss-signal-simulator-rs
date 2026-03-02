#!/usr/bin/env python3
"""
Верификация сгенерированного GNSS IF-файла.
Выполняет acquisition (захват спутников) для трёх систем:
  - GPS L1 C/A
  - BeiDou B1C (BOC(1,1))
  - Galileo E1 (BOC(1,1))

Если находим спутники с правильными PRN — сигнал сгенерирован корректно.
"""

import numpy as np
from scipy.fft import fft, ifft
import sys
import os
import re

# ============================================================
# Параметры сигнала
# ============================================================
SAMPLE_RATE = 5.0e6       # 5 MHz
IF_FREQ = 0.0             # Baseband (center = 1575.42 MHz)
CODE_RATE = 1.023e6       # GPS C/A chip rate
CODE_LENGTH = 1023        # GPS C/A chips per code period
CODE_PERIOD = 1e-3        # GPS C/A: 1 ms
SAMPLES_PER_MS = int(SAMPLE_RATE * CODE_PERIOD)  # 5000

# Acquisition parameters
DOPPLER_RANGE = 7000      # ±7 kHz search
DOPPLER_STEP = 250        # 250 Hz step
# Z-score thresholds: (peak - mean) / std — единый порог для всех систем
ZSCORE_THRESHOLD = 30     # Z-score порог (шум: ~5-10, сигнал: >50)
NON_COHERENT_GPS = 10     # non-coherent accumulations for GPS (1ms blocks)
NON_COHERENT_GAL = 20     # for Galileo (4ms blocks)
NON_COHERENT_BDS = 10     # for BeiDou (10ms blocks)

# ============================================================
# GPS C/A code generation (Gold codes)
# ============================================================

GPS_CA_TAPS = {
    1: (2, 6), 2: (3, 7), 3: (4, 8), 4: (5, 9), 5: (1, 9),
    6: (2, 10), 7: (1, 8), 8: (2, 9), 9: (3, 10), 10: (2, 3),
    11: (3, 4), 12: (5, 6), 13: (6, 7), 14: (7, 8), 15: (8, 9),
    16: (9, 10), 17: (1, 4), 18: (2, 5), 19: (3, 6), 20: (4, 7),
    21: (5, 8), 22: (6, 9), 23: (1, 3), 24: (4, 6), 25: (5, 7),
    26: (6, 8), 27: (7, 9), 28: (8, 10), 29: (1, 6), 30: (2, 7),
    31: (3, 8), 32: (4, 9),
}


def generate_ca_code(prn):
    """Generate GPS C/A code for given PRN (1-32). Returns ±1 array."""
    if prn not in GPS_CA_TAPS:
        return None
    g1 = np.ones(10, dtype=int)
    g2 = np.ones(10, dtype=int)
    tap1, tap2 = GPS_CA_TAPS[prn]
    code = np.zeros(CODE_LENGTH, dtype=np.int8)
    for i in range(CODE_LENGTH):
        g1_out = g1[9]
        g2_out = g2[tap1 - 1] ^ g2[tap2 - 1]
        code[i] = g1_out ^ g2_out
        g1_fb = g1[2] ^ g1[9]
        g1 = np.roll(g1, 1)
        g1[0] = g1_fb
        g2_fb = g2[1] ^ g2[2] ^ g2[5] ^ g2[7] ^ g2[8] ^ g2[9]
        g2 = np.roll(g2, 1)
        g2[0] = g2_fb
    return 1.0 - 2.0 * code.astype(np.float64)


# ============================================================
# BeiDou B1C code generation (Weil/Legendre sequences)
# ============================================================

# Parameters from prngenerate.rs
B1C_DATA_TRUNCATION = [
    699, 694, 7318, 2127, 715, 6682, 7850, 5495, 1162, 7682, 6792, 9973, 6596, 2092, 19, 10151,
    6297, 5766, 2359, 7136, 1706, 2128, 6827, 693, 9729, 1620, 6805, 534, 712, 1929, 5355, 6139,
    6339, 1470, 6867, 7851, 1162, 7659, 1156, 2672, 6043, 2862, 180, 2663, 6940, 1645, 1582, 951,
    6878, 7701, 1823, 2391, 2606, 822, 6403, 239, 442, 6769, 2560, 2502, 5072, 7268, 341,
]

B1C_DATA_PHASE_DIFF = [
    2678, 4802, 958, 859, 3843, 2232, 124, 4352, 1816, 1126, 1860, 4800, 2267, 424, 4192, 4333,
    2656, 4148, 243, 1330, 1593, 1470, 882, 3202, 5095, 2546, 1733, 4795, 4577, 1627, 3638, 2553,
    3646, 1087, 1843, 216, 2245, 726, 1966, 670, 4130, 53, 4830, 182, 2181, 2006, 1080, 2288, 2027,
    271, 915, 497, 139, 3693, 2054, 4342, 3342, 2592, 1007, 310, 4203, 455, 4318,
]


def generate_legendre_sequence(length):
    """Generate Legendre sequence of given length (matches Rust legendre_sequence)."""
    seq = np.zeros(length, dtype=np.int32)
    for i in range(1, length):
        seq[(i * i) % length] = 1
    return seq


def generate_b1c_weil(svid):
    """Generate BeiDou B1C data code for given SVID (1-63). Returns 0/1 array of length 10230."""
    if svid < 1 or svid > 63:
        return None
    trunc = B1C_DATA_TRUNCATION[svid - 1]
    phase = B1C_DATA_PHASE_DIFF[svid - 1]

    legendre = generate_legendre_sequence(10243)

    code = np.zeros(10230, dtype=np.int32)
    idx1 = trunc - 1
    idx2 = trunc + phase - 1

    for i in range(10230):
        if idx1 >= 10243:
            idx1 -= 10243
        if idx2 >= 10243:
            idx2 -= 10243
        code[i] = legendre[idx1] ^ legendre[idx2]
        idx1 += 1
        idx2 += 1

    return code


def apply_boc11(code_01):
    """Apply BOC(1,1) modulation to 0/1 code. Returns ±1 array of 2x length."""
    # Convert 0/1 to ±1: 0→+1, 1→-1 (matches Rust: prn!=0 → -1, else +1)
    code_pm = np.where(code_01 == 0, 1.0, -1.0)
    # BOC(1,1): each chip → [+chip, -chip] (subchip 1 is negated)
    boc = np.zeros(len(code_pm) * 2, dtype=np.float64)
    boc[0::2] = code_pm    # even subchips: same
    boc[1::2] = -code_pm   # odd subchips: negated
    return boc


def generate_b1c_boc(svid):
    """Generate BOC(1,1) modulated BeiDou B1C code. Returns ±1 array of 20460 subchips."""
    code = generate_b1c_weil(svid)
    if code is None:
        return None
    return apply_boc11(code)


# ============================================================
# Galileo E1 code generation (memory codes from file)
# ============================================================

_e1_memory_code = None


def load_e1_memory_codes():
    """Load Galileo E1 memory codes from Rust source file."""
    global _e1_memory_code
    if _e1_memory_code is not None:
        return _e1_memory_code

    filepath = os.path.join(os.path.dirname(__file__), 'src', 'memory_code_e1.rs')
    if not os.path.exists(filepath):
        print(f"  E1 memory code file not found: {filepath}")
        return None

    with open(filepath, 'r') as f:
        content = f.read()

    # Parse all hex values from the Rust array literal
    hex_values = re.findall(r'0x([0-9a-fA-F]+)', content)
    _e1_memory_code = np.array([int(h, 16) for h in hex_values], dtype=np.uint32)
    return _e1_memory_code


def extract_e1_code(svid, use_pilot=False):
    """Extract Galileo E1 code for given SVID (1-50). Returns 0/1 array of 4092 chips."""
    mem = load_e1_memory_codes()
    if mem is None:
        return None

    # Data: offset = (svid - 1) * 128, Pilot: offset = (svid + 49) * 128
    if use_pilot:
        offset = (svid + 49) * 128
    else:
        offset = (svid - 1) * 128

    if offset + 128 > len(mem):
        return None

    binary_seq = mem[offset:offset + 128]

    # Extract 4 sectors × 1023 chips (matches Rust get_memory_sequence)
    code = np.zeros(4092, dtype=np.int32)
    for sector in range(4):
        for j in range(1023):
            bit_index = j >> 5
            bit_offset = 31 - (j & 0x1f)
            word_index = sector * 32 + bit_index
            if word_index < len(binary_seq):
                code[sector * 1023 + j] = 1 if (int(binary_seq[word_index]) & (1 << bit_offset)) != 0 else 0

    return code


def generate_e1_boc(svid):
    """Generate BOC(1,1) modulated Galileo E1 data code. Returns ±1 array of 8184 subchips."""
    code = extract_e1_code(svid, use_pilot=False)
    if code is None:
        return None
    return apply_boc11(code)


# ============================================================
# Signal loading
# ============================================================

def load_iq8(filename, offset_ms=0, duration_ms=10):
    """Load IQ8 file. Returns complex64 numpy array."""
    n_samples = duration_ms * SAMPLES_PER_MS
    byte_offset = offset_ms * SAMPLES_PER_MS * 2

    filesize = os.path.getsize(filename)
    max_samples = filesize // 2
    if offset_ms * SAMPLES_PER_MS + n_samples > max_samples:
        n_samples = max_samples - offset_ms * SAMPLES_PER_MS

    raw = np.fromfile(filename, dtype=np.uint8, count=n_samples * 2, offset=byte_offset)
    iq = raw.astype(np.int8).astype(np.float64)
    I = iq[0::2]
    Q = iq[1::2]
    return I + 1j * Q


# ============================================================
# Resampling
# ============================================================

def resample_code(code, n_samples):
    """Resample code to n_samples points (nearest neighbor)."""
    indices = np.arange(n_samples) * len(code) / n_samples
    return code[indices.astype(int) % len(code)]


# ============================================================
# Generic FFT-based acquisition
# ============================================================

def acquire_satellite_generic(signal_blocks, local_code_resampled, doppler_range, doppler_step,
                              sample_rate, non_coherent_count):
    """
    Generic FFT-based parallel code phase search.
    signal_blocks: list of complex signal blocks (each = one code period)
    local_code_resampled: ±1 code resampled to match block length
    Returns (found, doppler_hz, code_phase, zscore).
    Uses z-score metric: (peak - mean) / std — independent of FFT size.
    """
    n = len(local_code_resampled)
    code_fft_conj = np.conj(fft(local_code_resampled))
    doppler_bins = np.arange(-doppler_range, doppler_range + 1, doppler_step)
    t = np.arange(n) / sample_rate

    best_zscore = 0
    best_doppler = 0
    best_code_phase = 0

    for doppler in doppler_bins:
        correlation_sum = np.zeros(n)
        count = min(non_coherent_count, len(signal_blocks))

        for blk in signal_blocks[:count]:
            carrier = np.exp(-1j * 2 * np.pi * (IF_FREQ + doppler) * t)
            baseband = blk[:n] * carrier
            sig_fft = fft(baseband)
            corr = np.abs(ifft(sig_fft * code_fft_conj)) ** 2
            correlation_sum += corr

        peak_idx = np.argmax(correlation_sum)
        peak_val = correlation_sum[peak_idx]

        # Z-score: исключаем область пика, считаем (peak - mean) / std
        mask = np.ones(n, dtype=bool)
        exclude = max(50, n // 100)
        mask[max(0, peak_idx - exclude):min(n, peak_idx + exclude)] = False
        noise = correlation_sum[mask]
        noise_mean = np.mean(noise)
        noise_std = np.std(noise)
        zscore = (peak_val - noise_mean) / noise_std if noise_std > 0 else 0

        if zscore > best_zscore:
            best_zscore = zscore
            best_doppler = doppler
            best_code_phase = peak_idx

    found = best_zscore > ZSCORE_THRESHOLD
    return found, best_doppler, best_code_phase, best_zscore


# ============================================================
# GPS L1 C/A Acquisition
# ============================================================

def acquire_gps(signal, prn):
    """Acquire GPS satellite. signal = complex array (at least 5ms)."""
    n = SAMPLES_PER_MS  # 5000 samples per 1ms
    ms_blocks = [signal[i * n:(i + 1) * n] for i in range(len(signal) // n)]
    ca_code = generate_ca_code(prn)
    if ca_code is None:
        return False, 0, 0, 0
    local_code = resample_code(ca_code, n)
    return acquire_satellite_generic(ms_blocks, local_code, DOPPLER_RANGE, DOPPLER_STEP,
                                     SAMPLE_RATE, NON_COHERENT_GPS)


# ============================================================
# BeiDou B1C Acquisition (BOC(1,1), 10ms code period)
# ============================================================

BDS_CODE_PERIOD_MS = 10
BDS_SAMPLES_PER_CODE = BDS_CODE_PERIOD_MS * SAMPLES_PER_MS  # 50000


def acquire_bds(signal, svid):
    """Acquire BeiDou B1C satellite. signal = complex array (at least 20ms)."""
    n = BDS_SAMPLES_PER_CODE  # 50000
    blocks = [signal[i * n:(i + 1) * n] for i in range(len(signal) // n)]
    boc_code = generate_b1c_boc(svid)
    if boc_code is None:
        return False, 0, 0, 0
    local_code = resample_code(boc_code, n)
    return acquire_satellite_generic(blocks, local_code, DOPPLER_RANGE, DOPPLER_STEP,
                                     SAMPLE_RATE, NON_COHERENT_BDS)


# ============================================================
# Galileo E1 Acquisition (BOC(1,1), 4ms code period)
# ============================================================

GAL_CODE_PERIOD_MS = 4
GAL_SAMPLES_PER_CODE = GAL_CODE_PERIOD_MS * SAMPLES_PER_MS  # 20000


def acquire_gal(signal, svid):
    """Acquire Galileo E1 satellite. signal = complex array (at least 12ms)."""
    n = GAL_SAMPLES_PER_CODE  # 20000
    blocks = [signal[i * n:(i + 1) * n] for i in range(len(signal) // n)]
    boc_code = generate_e1_boc(svid)
    if boc_code is None:
        return False, 0, 0, 0
    local_code = resample_code(boc_code, n)
    return acquire_satellite_generic(blocks, local_code, DOPPLER_RANGE, DOPPLER_STEP,
                                     SAMPLE_RATE, NON_COHERENT_GAL)


# ============================================================
# Signal quality analysis
# ============================================================

def analyze_signal_stats(signal):
    """Basic signal statistics."""
    I = np.real(signal)
    Q = np.imag(signal)

    print("=" * 60)
    print("  СТАТИСТИКА СИГНАЛА")
    print("=" * 60)
    print(f"  Сэмплов:     {len(signal):,}")
    print(f"  Длительность: {len(signal)/SAMPLE_RATE:.3f} с")
    print(f"  I: mean={np.mean(I):.3f}, std={np.std(I):.3f}, "
          f"min={np.min(I):.0f}, max={np.max(I):.0f}")
    print(f"  Q: mean={np.mean(Q):.3f}, std={np.std(Q):.3f}, "
          f"min={np.min(Q):.0f}, max={np.max(Q):.0f}")
    print(f"  Мощность I:  {np.mean(I**2):.2f}")
    print(f"  Мощность Q:  {np.mean(Q**2):.2f}")
    print(f"  Полная мощн.: {np.mean(np.abs(signal)**2):.2f}")
    print()


def compute_spectrum(signal, n_fft=4096):
    """Compute power spectrum in dB."""
    n_avg = min(len(signal) // n_fft, 50)
    psd = np.zeros(n_fft)
    for i in range(n_avg):
        seg = signal[i * n_fft:(i + 1) * n_fft]
        window = np.hanning(n_fft)
        psd += np.abs(fft(seg * window)) ** 2
    psd /= n_avg
    psd_db = 10 * np.log10(psd + 1e-20)
    freqs = np.fft.fftfreq(n_fft, 1.0 / SAMPLE_RATE) / 1e6
    return np.fft.fftshift(freqs), np.fft.fftshift(psd_db)


# ============================================================
# Run acquisition for one system
# ============================================================

def run_acquisition(system_name, signal, svid_range, acquire_fn):
    """Run acquisition for a GNSS system and return found satellites."""
    print(f"\n{'=' * 60}")
    print(f"  {system_name} ACQUISITION")
    print(f"{'=' * 60}")
    print(f"  Допплер: ±{DOPPLER_RANGE} Гц, шаг {DOPPLER_STEP} Гц")
    print(f"  Порог Z-score: {ZSCORE_THRESHOLD}")
    print()

    found_sats = []
    all_results = []
    total = len(svid_range)
    for idx, svid in enumerate(svid_range):
        print(f"\r  Поиск {system_name} SV{svid:02d}... ({idx+1}/{total})", end='', flush=True)
        found, doppler, code_phase, zscore = acquire_fn(signal, svid)
        all_results.append((svid, doppler, code_phase, zscore))
        if found:
            found_sats.append((svid, doppler, code_phase, zscore))

    print(f"\r  {'':60}")  # clear line

    # Показываем top-10 по z-score для диагностики
    top10 = sorted(all_results, key=lambda x: -x[3])[:10]
    print(f"  {'SV':>4}  {'Допплер':>10}  {'Z-score':>10}  {'Статус'}")
    print(f"  {'─'*4}  {'─'*10}  {'─'*10}  {'─'*12}")
    for svid, doppler, cp, zs in top10:
        status = "НАЙДЕН" if zs > ZSCORE_THRESHOLD else "шум"
        print(f"  {svid:4d}  {doppler:>+8.0f} Гц  {zs:>9.1f}  {status}")

    print(f"\n  РЕЗУЛЬТАТ: Найдено {len(found_sats)} {system_name} спутников")
    return found_sats


# ============================================================
# Main
# ============================================================

def main():
    filename = "generated_files/GPS_BDS_GAL_triple_10s.C8"
    if len(sys.argv) > 1:
        filename = sys.argv[1]

    if not os.path.exists(filename):
        print(f"Файл не найден: {filename}")
        return

    filesize = os.path.getsize(filename)
    total_ms = filesize // (SAMPLES_PER_MS * 2)
    print(f"\nФайл: {filename}")
    print(f"Размер: {filesize / 1e6:.1f} МБ")
    print(f"Длительность: {total_ms} мс ({total_ms/1000:.1f} с)")
    print(f"Sample rate: {SAMPLE_RATE/1e6:.1f} МГц")

    # Load enough data for non-coherent accumulation
    # BDS: 10 * 10ms = 100ms, GAL: 20 * 4ms = 80ms, GPS: 10 * 1ms = 10ms
    load_ms = min(500, total_ms)
    print(f"\nЗагружаю первые {load_ms} мс для анализа...")
    signal = load_iq8(filename, offset_ms=0, duration_ms=load_ms)

    # Signal statistics
    analyze_signal_stats(signal)

    # Spectrum
    freqs, psd_db = compute_spectrum(signal)

    # === GPS L1 C/A ===
    gps_sats = run_acquisition("GPS L1CA", signal, range(1, 33), acquire_gps)

    # === BeiDou B1C ===
    bds_sats = run_acquisition("BeiDou B1C", signal, range(1, 64), acquire_bds)

    # === Galileo E1 ===
    gal_sats = run_acquisition("Galileo E1", signal, range(1, 37), acquire_gal)

    # === Summary ===
    print(f"\n{'=' * 60}")
    print(f"  ИТОГОВАЯ СВОДКА")
    print(f"{'=' * 60}")
    print(f"  GPS L1CA:    {len(gps_sats)} спутников")
    print(f"  BeiDou B1C:  {len(bds_sats)} спутников")
    print(f"  Galileo E1:  {len(gal_sats)} спутников")
    print(f"  ВСЕГО:       {len(gps_sats) + len(bds_sats) + len(gal_sats)} спутников")

    if gps_sats:
        print(f"\n  GPS: {', '.join(f'G{s[0]:02d}(z={s[3]:.0f})' for s in sorted(gps_sats, key=lambda x: -x[3]))}")
    if bds_sats:
        print(f"  BDS: {', '.join(f'C{s[0]:02d}(z={s[3]:.0f})' for s in sorted(bds_sats, key=lambda x: -x[3]))}")
    if gal_sats:
        print(f"  GAL: {', '.join(f'E{s[0]:02d}(z={s[3]:.0f})' for s in sorted(gal_sats, key=lambda x: -x[3]))}")

    # === Plot ===
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        fig.suptitle(f'GNSS Signal Verification: {os.path.basename(filename)}', fontsize=14)

        # 1. Spectrum
        ax = axes[0, 0]
        ax.plot(freqs, psd_db, linewidth=0.5)
        ax.set_xlabel('Frequency offset (MHz)')
        ax.set_ylabel('Power (dB)')
        ax.set_title('Power Spectrum')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-2.5, 2.5)

        # 2. I/Q histogram
        ax = axes[0, 1]
        I_vals = np.real(signal[:50000]).astype(int)
        Q_vals = np.imag(signal[:50000]).astype(int)
        bins = np.arange(-130, 130, 1)
        ax.hist(I_vals, bins=bins, alpha=0.6, label='I', density=True)
        ax.hist(Q_vals, bins=bins, alpha=0.6, label='Q', density=True)
        ax.set_xlabel('Sample value')
        ax.set_ylabel('Density')
        ax.set_title('I/Q Histogram')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # 3. Constellation diagram
        ax = axes[0, 2]
        n_plot = min(5000, len(signal))
        ax.scatter(np.real(signal[:n_plot]), np.imag(signal[:n_plot]), s=0.5, alpha=0.3)
        ax.set_xlabel('I')
        ax.set_ylabel('Q')
        ax.set_title('I/Q Constellation')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)

        # 4. GPS correlation (strongest satellite)
        ax = axes[1, 0]
        if gps_sats:
            best = sorted(gps_sats, key=lambda x: -x[3])[0]
            ca_code = generate_ca_code(best[0])
            local_code = resample_code(ca_code, SAMPLES_PER_MS)
            code_fft_conj = np.conj(fft(local_code))
            t = np.arange(SAMPLES_PER_MS) / SAMPLE_RATE
            carrier = np.exp(-1j * 2 * np.pi * best[1] * t)
            corr_sum = np.zeros(SAMPLES_PER_MS)
            for i in range(min(NON_COHERENT_GPS, load_ms)):
                blk = signal[i * SAMPLES_PER_MS:(i + 1) * SAMPLES_PER_MS]
                baseband = blk * carrier
                corr = np.abs(ifft(fft(baseband) * code_fft_conj)) ** 2
                corr_sum += corr
            code_axis = np.arange(SAMPLES_PER_MS) * CODE_LENGTH / SAMPLES_PER_MS
            ax.plot(code_axis, corr_sum / np.max(corr_sum), linewidth=0.5)
            ax.set_title(f'GPS PRN {best[0]} (Doppler {best[1]:+.0f} Hz, SNR {best[3]:.1f})')
        else:
            ax.text(0.5, 0.5, 'No GPS satellites', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('GPS Correlation')
        ax.set_xlabel('Code phase (chips)')
        ax.set_ylabel('Normalized correlation')
        ax.grid(True, alpha=0.3)

        # 5. BeiDou correlation (strongest satellite)
        ax = axes[1, 1]
        if bds_sats:
            best = sorted(bds_sats, key=lambda x: -x[3])[0]
            boc_code = generate_b1c_boc(best[0])
            local_code = resample_code(boc_code, BDS_SAMPLES_PER_CODE)
            code_fft_conj = np.conj(fft(local_code))
            t = np.arange(BDS_SAMPLES_PER_CODE) / SAMPLE_RATE
            carrier = np.exp(-1j * 2 * np.pi * best[1] * t)
            corr_sum = np.zeros(BDS_SAMPLES_PER_CODE)
            n_blk = len(signal) // BDS_SAMPLES_PER_CODE
            for i in range(min(NON_COHERENT_BDS, n_blk)):
                blk = signal[i * BDS_SAMPLES_PER_CODE:(i + 1) * BDS_SAMPLES_PER_CODE]
                baseband = blk * carrier
                corr = np.abs(ifft(fft(baseband) * code_fft_conj)) ** 2
                corr_sum += corr
            code_axis = np.arange(BDS_SAMPLES_PER_CODE) * 10230 / BDS_SAMPLES_PER_CODE
            ax.plot(code_axis, corr_sum / np.max(corr_sum), linewidth=0.3)
            ax.set_title(f'BDS C{best[0]:02d} (Doppler {best[1]:+.0f} Hz, SNR {best[3]:.1f})')
        else:
            ax.text(0.5, 0.5, 'No BeiDou satellites', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('BeiDou Correlation')
        ax.set_xlabel('Code phase (chips)')
        ax.set_ylabel('Normalized correlation')
        ax.grid(True, alpha=0.3)

        # 6. Galileo correlation (strongest satellite)
        ax = axes[1, 2]
        if gal_sats:
            best = sorted(gal_sats, key=lambda x: -x[3])[0]
            boc_code = generate_e1_boc(best[0])
            local_code = resample_code(boc_code, GAL_SAMPLES_PER_CODE)
            code_fft_conj = np.conj(fft(local_code))
            t = np.arange(GAL_SAMPLES_PER_CODE) / SAMPLE_RATE
            carrier = np.exp(-1j * 2 * np.pi * best[1] * t)
            corr_sum = np.zeros(GAL_SAMPLES_PER_CODE)
            n_blk = len(signal) // GAL_SAMPLES_PER_CODE
            for i in range(min(NON_COHERENT_GAL, n_blk)):
                blk = signal[i * GAL_SAMPLES_PER_CODE:(i + 1) * GAL_SAMPLES_PER_CODE]
                baseband = blk * carrier
                corr = np.abs(ifft(fft(baseband) * code_fft_conj)) ** 2
                corr_sum += corr
            code_axis = np.arange(GAL_SAMPLES_PER_CODE) * 4092 / GAL_SAMPLES_PER_CODE
            ax.plot(code_axis, corr_sum / np.max(corr_sum), linewidth=0.3)
            ax.set_title(f'GAL E{best[0]:02d} (Doppler {best[1]:+.0f} Hz, SNR {best[3]:.1f})')
        else:
            ax.text(0.5, 0.5, 'No Galileo satellites', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Galileo Correlation')
        ax.set_xlabel('Code phase (chips)')
        ax.set_ylabel('Normalized correlation')
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plot_path = 'generated_files/signal_verification.png'
        plt.savefig(plot_path, dpi=150)
        print(f"\n  График сохранён: {plot_path}")
        plt.close()

    except Exception as e:
        print(f"\n  Не удалось сохранить график: {e}")

    print("\nГотово!")


if __name__ == '__main__':
    main()
