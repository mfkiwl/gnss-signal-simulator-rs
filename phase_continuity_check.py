#!/usr/bin/env python3
"""Carrier-phase continuity checker for GPS L1CA in an IQ8 file.

Acquires the strongest visible GPS L1CA satellite, then tracks its carrier phase
millisecond-by-millisecond by removing a *continuous* carrier replica (IF + Doppler
referenced to absolute sample time). Any residual phase step at the 1 ms boundaries is
a generator-side discontinuity. A clean generator gives ~0 step after detrending; the
audit H1 bug (lost IF fractional cycles) shows a ~0.5 cycle/ms sawtooth.

Usage: python phase_continuity_check.py <iq.C8> --preset presets/gps_gal_l1.json [--prn N]
"""
import argparse
import json
import sys
import numpy as np

# --- GPS L1CA code generation (IS-GPS-200) ---
G2_DELAY = {
    1: 5, 2: 6, 3: 7, 4: 8, 5: 17, 6: 18, 7: 139, 8: 140, 9: 141, 10: 251,
    11: 252, 12: 254, 13: 255, 14: 256, 15: 257, 16: 258, 17: 469, 18: 470,
    19: 471, 20: 472, 21: 473, 22: 474, 23: 509, 24: 512, 25: 513, 26: 514,
    27: 515, 28: 516, 29: 859, 30: 860, 31: 861, 32: 862,
}


def l1ca_code(prn):
    g1 = [1] * 10
    g2 = [1] * 10
    out = np.empty(1023, dtype=np.int8)
    delay = G2_DELAY[prn]
    g2_hist = []
    for i in range(1023 + delay):
        g1_out = g1[9]
        g2_out = g2[9]
        g2_hist.append(g2_out)
        nb1 = g1[2] ^ g1[9]
        nb2 = g2[1] ^ g2[2] ^ g2[5] ^ g2[7] ^ g2[8] ^ g2[9]
        g1 = [nb1] + g1[:9]
        g2 = [nb2] + g2[:9]
        if i >= delay:
            ca = g1_out ^ g2_hist[i - delay]
            out[i - delay] = 1 - 2 * ca  # 0/1 -> +1/-1
    return out.astype(np.float64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("iq")
    ap.add_argument("--preset", required=True)
    ap.add_argument("--prn", type=int, default=0)
    ap.add_argument("--ms", type=int, default=400, help="ms to track")
    args = ap.parse_args()

    preset = json.load(open(args.preset))
    out = preset["output"]
    fs = float(out["sampleFreq"]) * 1e6
    center = float(out["centerFreq"]) * 1e6
    if_freq = (1575.42e6 - center)  # GPS L1 IF offset (Hz)
    spms = int(round(fs / 1000.0))  # samples per ms

    raw = np.fromfile(args.iq, dtype=np.int8)
    iq = raw[0::2].astype(np.float64) + 1j * raw[1::2].astype(np.float64)
    total_ms = len(iq) // spms
    use_ms = min(args.ms, total_ms)

    # upsampled code helper
    def code_samples(code, n, fs, chip_rate=1.023e6):
        idx = (np.arange(n) * (chip_rate / fs)).astype(np.int64) % 1023
        return code[idx]

    prns = [args.prn] if args.prn else range(1, 33)
    best = None
    acq_ms = min(20, total_ms)  # non-coherent integration length
    tn = np.arange(spms) / fs
    for prn in prns:
        code = l1ca_code(prn)
        c = code_samples(code, spms, fs)
        cf = np.conj(np.fft.fft(c))
        for dopp in np.arange(-5000, 5001, 250):
            f = if_freq + dopp
            acc = np.zeros(spms)
            for kb in range(acq_ms):
                n0 = kb * spms
                s = iq[n0:n0 + spms] * np.exp(-1j * 2 * np.pi * f * (n0 / fs + tn))
                corr = np.fft.ifft(np.fft.fft(s) * cf)
                acc += np.abs(corr) ** 2  # non-coherent power sum
            mag = np.sqrt(acc)
            z = mag.max() / (mag.mean() + 1e-12)
            if best is None or z > best[0]:
                best = (z, prn, dopp, int(np.argmax(mag)), f)

    z, prn, dopp, code_phase, f = best
    print(f"Acquired GPS L1CA PRN{prn:02d}: Doppler={dopp:+.0f} Hz, IF={if_freq:+.0f} Hz, "
          f"carrier={f:+.0f} Hz, code_phase={code_phase}, z={z:.1f}")
    if z < 8:
        print("  (weak acquisition — phase result may be unreliable)")

    # --- track per-ms carrier phase with a CONTINUOUS replica (absolute time) ---
    code = l1ca_code(prn)
    phases = []
    for k in range(use_ms):
        n0 = k * spms
        n = np.arange(n0, n0 + spms)
        t = n / fs
        carrier = np.exp(-1j * 2 * np.pi * f * t)  # continuous across ms boundaries
        c = code_samples(code, spms, fs)
        # align code phase
        c = np.roll(c, code_phase)
        corr = np.sum(iq[n0:n0 + spms] * carrier * c)
        phases.append(np.angle(corr) / (2 * np.pi))
    phases = np.array(phases)

    # unwrap, remove linear trend (residual freq error), measure boundary steps
    unwrapped = np.unwrap(phases * 2 * np.pi) / (2 * np.pi)
    k = np.arange(len(unwrapped))
    slope, intercept = np.polyfit(k, unwrapped, 1)
    residual = unwrapped - (slope * k + intercept)
    steps = np.abs(np.diff(residual))
    print(f"  tracked {use_ms} ms | residual freq slope={slope*1000:+.1f} Hz")
    print(f"  max per-ms phase step (detrended) = {steps.max():.4f} cycles")
    print(f"  mean per-ms phase step            = {steps.mean():.4f} cycles")
    print(f"  std of detrended phase            = {residual.std():.4f} cycles")
    # half-cycle sawtooth signature (H1): alternating ~0.5 steps
    near_half = np.mean((steps > 0.3) & (steps < 0.7))
    print(f"  fraction of ~0.5-cycle steps      = {near_half:.2%}")


if __name__ == "__main__":
    main()
