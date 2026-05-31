# Полный аудит кода gnss_rust — 2026-05-31

Многоагентный аудит: 18 finder-агентов по модулям/кросс-срезам, каждая High/Medium/Critical находка прошла состязательную верификацию (скептик пытался опровергнуть). 90 агентов, 106 находок, 96 подтверждено, 8 отклонено.

**Severity:** 1 Critical · 23 High · 31 Medium · 44 Low


---

# Аудит корректности GNSS IF-симулятора (gnss_rust)

Финальный отчёт ведущего ревьюера. Все находки подтверждены чтением исходного кода, сверкой с GNSS ICD и (где применимо) численным воспроизведением. Сгруппировано по серьёзности, дубликаты с общей первопричиной объединены.

---

## Critical

### C1. F/NAV свёрточный кодер: неверный полином G1 и отсутствует инверсия G2
**Файл:** `src/fnavbit.rs:94-102`

F/NAV FEC использует `count1(conv_state & 0x75)` для G1 и `count1(conv_state & 0x5B)` для G2 без инверсии. Galileo K=7 r=1/2 код (общий для F/NAV и I-NAV) требует G1=171oct, G2=133oct с обязательной инверсией выхода G2. Эталонный I-NAV кодер в этом же репозитории (`src/inavbit.rs:980-995`, инверсия `^0x5` на строке 942) реализует правильный код. При данном порядке регистра (`conv_state = (conv_state<<1)|bit`, новейший бит в LSB) корректные маски — **G1=0x4F**, **G2=0x6D с инверсией**, что воспроизводит ICD-кодер 200/200 символов. Текущие маски `0x75`/`0x5B` расходятся со стандартом примерно на половине символов; FNavBit — живой генератор nav-данных для Galileo E5a (`nav_data.rs:351`, диспетчеризация через `satellite_signal.rs:391`). Любой приёмник, выполняющий Viterbi-декодирование E5a, не сможет восстановить навигационное сообщение.

**Важно:** одноточечный фикс из исходной заявки (`& 0x5B ^ 1` при сохранении G1=0x75) НЕ исправляет кодер — нужно заменить ОБЕ маски.

**Фикс:**
```rust
// G1 = 171 oct, G2 = 133 oct (инвертирован) — как в inavbit.rs
encoded_symbols[i * 2]     = count1(conv_state & 0x4F);
encoded_symbols[i * 2 + 1] = count1(conv_state & 0x6D) ^ 1;
```
Исправить также вводящий в заблуждение комментарий на строке 94.

---

## High

### H1. Целый класс багов фазы несущей: дробные циклы IF теряются на каждой границе мс
**Файлы:** `src/sat_if_signal.rs:654-669` (и AVX-путь `494-501`)

*Объединяет два дубликата заявки — одна первопричина.*

В `get_if_sample_cached` фаза за один мс продвигается на `doppler_cycles_per_ms + if_freq/1000` циклов (строка 656), но межмиллисекундный якорь обновляется только доплеровской частью: `self.start_carrier_phase -= doppler_cycles_per_ms;` (строка 669). Вклад IF `frac(if_freq/1000)` никогда не переносится в стартовую фазу следующего мс, что вносит детерминированный скачок фазы `frac(if_freq/1000)` циклов каждую миллисекунду. Численно подтверждено: при `if_freq` кратном 1000 Гц скачок 0.000; при `if_freq=1234567` Гц — 0.433 цикла/мс. Для пресета `presets/quad_l1g1.json` (centerFreq=1582.2105 МГц) GPS/BDS/GAL получают `if_freq=-6790500` Гц → скачок **0.5 цикла (180°) каждую мс**; GLONASS FDMA (k=0) при том же пресете даёт `if_freq=19789500` Гц → также 0.5 цикла/мс. Компенсация GLONASS на строках 661-666 закрывает только случай нечётных k и ровно 0.5; GPS/BDS/GAL и нецелые дроби не покрываются. Эталон `src/gps_pilot/mod.rs:286` аккумулирует ПОЛНУЮ частоту непрерывно (`carrier_phase += freq*dt`) и документирует 0.000000 цикла макс. скачка. Это разрушает непрерывность фазы / PLL-tracking; соответствует документированному результату «GLONASS только 5/7 детектировано».

**Фикс:** аккумулировать полное продвижение в якорь — `self.start_carrier_phase -= doppler_cycles_per_ms + self.if_freq as f64/1000.0;` — или перейти на непрерывную посэмпловую аккумуляцию из `gps_pilot/mod.rs`, после чего GLONASS-спецслучай half-cycle можно удалить.

---

### H2. GLONASS G3: IF-частота портится до −CenterFreq (G3 генерируется на ~1.2 ГГц мимо baseband)
**Файл:** `src/sat_if_signal.rs:360-364`

`get_signal_center_freq()` в ветке GLONASS обрабатывает только `SIGNAL_INDEX_G1` и `SIGNAL_INDEX_G2`; `SIGNAL_INDEX_G3` (=26) попадает в `_ => 0.0`. Этот 0.0 используется в `update_satellite_params` (315-316) как `if_freq = (signal_center_freq - output_center_freq) as i32`, перезаписывая корректное значение от `create_glonass_signals`. При пресете `glo_g3.json` (centerFreq 1202.025 МГц) получается `(0.0 - 1202025000.0) = -1202025000`, и несущая каждого G3-сэмпла оказывается на ~1.2 ГГц мимо baseband — вне 12-МГц полосы, G3 недетектируем. `PARAM_UPDATE_INTERVAL_MS=1` срабатывает на первом же мс. То, что это баг, а не намеренное поведение, доказывает `src/satellite_param.rs:494-501` (`get_wave_length`), где идентичный маппинг явно содержит `SIGNAL_INDEX_G3 => FREQ_GLO_G3, // G3 uses CDMA, not FDMA`.

**Фикс:** добавить в ветку GLONASS `SIGNAL_INDEX_G3 => FREQ_GLO_G3,` (G3/L3OC — CDMA, без FDMA-сдвига).

---

### H3. GLONASS L2P: code phase уничтожается каждую мс из-за `rem_euclid(pilot_length=1)`
**Файл:** `src/sat_if_signal.rs:786-787` (с `224-227`)

Для GPS L2P `new()` ставит `data_length=20460`, но принудительно `pilot_length=1` (строка 226 `pl = 1`). В конце каждого мс трекер обёртывается `start_code_phase.rem_euclid(1.0)`, схлопывая счётчик в `[0,1)`. `base_chip_offset = start_code_phase`, поэтому за каждый мс индекс чипа проходит лишь chips 0..~10230 из 20460-чиповой эпохи и рестартует у нуля. Вторая половина каждой L2P-эпохи никогда не излучается, код не пересекает 2-мс границу (20460 чипов). `init_state` корректно якорит первую мс через реальный `pilot_period=2`, но `update_satellite_params` намеренно не реанкорит, так что с мс 1 всё ломается. AVX-путь редиректит всё с `data_length != 1023` сюда, поэтому L2P всегда идёт через баг. `presets/gps_l2p.json` запускает этот путь — IQ не является валидным P-кодом.

**Фикс:** обёртывать по реальному периоду: `let wrap_len = if self.pilot_length > 1 { self.pilot_length } else { self.data_length };` и `rem_euclid(wrap_len as f64)`; либо в L2P-спецслучае ставить `pl = 20460`.

---

### H4. GLONASS-пропагатор: двойное масштабирование км→м (лишний ×1000) в `glonass_sat_pos_speed_eph`
**Файл:** `src/satellite_param.rs:849-859` (используется `get_glonass_satellite_param`, вызовы `ifdatagen.rs:3373,3768`)

Локальный RK4-пропагатор строит начальное состояние как `eph.x * 1000.0` (и так же для y/z/v/a), хотя живой парсер `parse_glonass_ephemeris_correct` (`json_interpreter.rs:3382`) уже конвертирует в метры (`x: data[3] * 1e3`). RK4-derivatives используют SI-константы (`PZ90_GM=3.986004418e14`, `PZ90_AE=6378136.0`), которым нужна позиция в метрах. Раздутие состояния до ~2.6e10 делает гравитационное ускорение `-GM/r³` слабее на ~10⁶, интегрирование вырождается в прямую. Численно: ошибка позиции 1.1 км при dt=60с, 27.6 км при dt=300с, 248 км при dt=900с. Поскольку tb GLONASS привязан к 15-мин окну, delta_t регулярно достигает минут → неверные Elevation/Azimuth/TravelTime/Doppler. Эталон `coordinate.rs:451-461` использует `eph.x` напрямую без перемасштабирования и применяет вращение скорости в инерциальную систему (`eph.vx - PZ90_OMEGDOTE*eph.y`), которое локальная версия теряет. Дубликат той же ошибки — `ifdatagen.rs:4785`.

**Фикс:** убрать `* 1000.0` на входе и `/ 1000.0` на выходе; либо делегировать GLONASS в `crate::coordinate::glonass_sat_pos_speed_eph`.

---

### H5. GLONASS-видимость использует no-op заглушку позиции, игнорирующую время симуляции
**Файл:** `src/ifdatagen.rs:4775-4792`

`IFDataGen::glonass_sat_pos_speed_eph` (вызов из цикла видимости `1991`) полностью игнорирует аргумент `_transmit_time` и возвращает сырую референсную позицию эфемериды (km→m) без какой-либо пропагации орбиты (комментарий: «In a real implementation, this would use Runge-Kutta integration»). GLONASS-спутники выбираются как видимые по позиции на эпохе tb, а не на эпохе симуляции; при spacing tb ~30 мин |t−tb| доходит до ~900с, ~3.9 км/с → тысячи км немоделируемого движения, изменение elevation на несколько градусов у границы маски. При этом сигнал генерируется через настоящий RK4 (`get_glonass_satellite_param`, `satellite_param.rs:833`), так что видимый и генерируемый наборы считаются по несогласованным позициям. GPS/BeiDou/Galileo-видимость (циклы `1757/1860/1930`) корректно используют пропагирующий `gps_sat_pos_speed_eph`.

**Фикс:** заменить заглушку на `crate::coordinate::glonass_sat_pos_speed_eph(transmit_time, ...)`.

---

### H6. `utc_to_gps_time` дважды учитывает дробные секунды
**Файл:** `src/gnsstime.rs:236-243`

`utc_to_gps_time` пишет `MilliSeconds = seconds_in_week*1000 + ((Second%1.0)*1000.0) as u32` (строка 237) И одновременно `SubMilliSeconds = Second%1.0` (строка 243). Контракт кодбейса (`json_interpreter.rs:547-549`, `satellite_param.rs:257` где `time = (MilliSeconds + SubMilliSeconds)/1000.0`): `MilliSeconds` — целые мс, `SubMilliSeconds` — остаток в долях ОДНОЙ мс `[0,1)`. Здесь полусекунда попадает и в `MilliSeconds` (+500), и повторно как +0.5 мс. `get_transmit_time` заимствует лишь одну мс и не может размотать это; `add_milliseconds` пропускает `SubMilliSeconds` без нормализации. Для `Second=N.5` получается спурьёзный ~0.5 мс сдвиг на мастер-базе времени (~150 км ошибки псевдодальности). Латентно, т.к. все пресеты используют целые секунды; `test_utc_to_gps_conversion` использует `Second=0.0`.

**Фикс:** зеркалить `json_interpreter.rs`: `let ms_f = Second*1000.0; let ms_int = ms_f as i32; MilliSeconds = seconds_in_week*1000 + ms_int; SubMilliSeconds = ms_f - ms_int as f64;`. Добавить тест с дробной стартовой секундой.

---

### H7. RINEX: off-by-one в защите полей фиксированной ширины (строгий `>` теряет последнее поле обрезанных строк)
**Файлы:** `src/json_interpreter.rs:3277,3289,3296` (`read_contents_data`); `3202,3209` (`read_contents_time`)

19-символьные поля заканчиваются на колонках 23/42/61/80. Защиты используют строгий `>` (`if line.len() > 23 {...}`). Когда RINEX-писатель обрезает trailing whitespace, строка, кончающаяся ровно на границе поля, имеет длину == границе, и `boundary > boundary` ложно — поле молча остаётся 0.0. Не латентно: репозиторная фикстура `BRDC00IGS_R_20251560000_01D_MN.rnx` (используется `cargo run --bin nav_diagnostics`) содержит 11235 строк длиной 23, 1568 длиной 42, 7004 длиной 61. Конкретно orbit-5 у Galileo E02 имеет длину 61, поле data[21] (GST week) теряется → `eph.week=0` → `corrected_week=1024`, эфемерида улетает на ~45 лет и отбрасывается 3-часовым фильтром (правдоподобно объясняет недосчёт Galileo). Поле 4 (граница 80) восстанавливается fallback-веткой `else if line.len() > 61`, так что реальные потери — поля 1-3 (week, IDOT/SISA, transmission time).

**Фикс:** заменить на `>=` (`>= 23`, `>= 42`, `>= 61`) и клампить верхний индекс среза: `&line[4..23.min(line.len())]`.

---

### H8. BeiDou: TGD2 теряется, IODC портится (GPS-парсер читает data[26] как IODC, tgd2 не присваивается)
**Файл:** `src/json_interpreter.rs:2671-2675, 3001, 3138`

`parse_beidou_ephemeris` делегирует в `parse_gps_ephemeris`, затем копирует `tgd2: gps_eph.tgd2`. Но `eph.tgd2 =` присваивается ТОЛЬКО в Galileo-парсере (строка 3138), никогда в GPS-парсере, поэтому `tgd2=0.0` (Default) для каждой BeiDou-записи. Для BeiDou orbit-6 — `[SV accuracy, SatH1, TGD1, TGD2]`, т.е. data[26]=TGD2; GPS-парсер делает `eph.iodc = data[26] as u16`, загружая TGD2 (малое отрицательное, напр. -9.7e-9) в iodc, при `f64→u16 as` (насыщающий) → 0. Проверено по `rinex_v3_20251560000.rnx` C01: data[26]=-9.7e-9. `d1d2navbit.rs:351,506` кодирует `tgd2 * 1e10` в D1/D2 nav-сообщение, поэтому излучаемый сигнал несёт TGD2=0 вместо реальной групповой задержки B2I.

**Фикс:** в ветке BeiDou явно `tgd1 = data[25]`, `tgd2 = data[26]`, и брать iodc/aodc из настоящего поля AODC (orbit-7 / data[28]) вместо data[26].

---

### H9. L5 CNAV MT30: неверные масштабы часовых членов (af0 2⁻³⁴, af1 2⁻⁴⁶, af2 2⁻⁵⁹)
**Файл:** `src/l5cnavbit.rs:393-402`

`compose_eph_words_static` кодирует MT30 с `af0*2³⁴`, `af1*2⁴⁶`, `af2*2⁵⁹`. IS-GPS-705 MT30 (идентично IS-GPS-200): af0 — 26 бит scale 2⁻³⁵, af1 — 20 бит 2⁻⁴⁸, af2 — 10 бит 2⁻⁶⁰. Эталон `cnavbit.rs:394-400` использует `-35/-48/-60`. L5 ошибается ровно в 2×/4×/2×. Приёмник, декодирующий L5 CNAV, получает clock bias/drift/drift-rate с неверным масштабом (af0 — доминирующий метровый член).

**Фикс:** `2.0f64.powi(35)` для af0, `powi(48)` для af1, `powi(60)` для af2 — как в `cnavbit.rs`.

---

### H10. L5 CNAV MT10/MT11: кеплеровы углы 2⁻³¹/32 бит вместо 2⁻³²/33 бит; кодируется sqrtA вместо ΔA
**Файл:** `src/l5cnavbit.rs:344-355, 362-372`

M0, omega0, i0, w кодируются как 32-битные при scale 2⁻³¹ семициклов; sqrtA пакуется при 2⁻¹⁹. CNAV MT10/MT11 (IS-GPS-705/200) задаёт эти углы как 33-битные при 2⁻³² семициклов (эталон `cnavbit.rs:309,329,344,354` через `unscale_long(.../PI,-32)` с явным 33-м знаковым битом). L5 теряет MSB и половинит разрешение. Кроме того, CNAV MT10 передаёт ΔA = axis − A_REF (26 бит, 2⁻⁹ м) + axis_dot, а не sqrtA. Поля `axis`/`axis_dot` уже есть в `GpsEphemeris` и используются рабочим L2C-путём. Структура полей L5 MT10/MT11 также неполна (нет Health, top, URA, axis_dot, delta_n). Эфемерида декодируется в неверные орбитальные параметры.

**Фикс:** как `cnavbit.rs` — `unscale_long(.../PI,-32)` с сохранением 33-го бита; кодировать `ΔA = axis - A_REF` при 2⁻⁹ в MT10 с axis_dot.

---

### H11. L5 CNAV MT30: TGD/ISC кодируются как 10 бит 2⁻³² вместо 13 бит 2⁻³⁵
**Файл:** `src/l5cnavbit.rs:414-424`

TGD и ISC-члены пакуются `tgd * 2³²`, маска `& 0x3FF` (10 бит). IS-GPS-705 MT30: TGD/ISC — 13-битные two's-complement при LSB 2⁻³⁵ с (эталон `cnavbit.rs:404-415` через `unscale_int(...,-35)` в 13-битные поля). Масштаб 2⁻³² против 2⁻³⁵ — 8× грубее, а 10 бит против 13 размещают биты в неверных позициях кадра. Весь MT30-блок L5 дефектен (см. также H9).

**Фикс:** кодировать TGD и ISC при 2⁻³⁵ в 13-битных полях, как `cnavbit.rs`.

---

### H12. CNAV MT33 (UTC): data[7] дублирует utc_message[2] вместо [3], теряя DN/TLSF/WNLSF-LSB
**Файл:** `src/cnavbit.rs:541-543`

Для MT33 кодер пишет `data[6]=utc_message[2]; data[7]=utc_message[2]; data[8]=0;` — дублирует слово 2 и никогда не копирует `utc_message[3]`. `set_iono_utc` (268-272) пакует UTC в четыре слова: `utc_message[3]` несёт WNLSF-LSB (бит31), DN (биты27-30), TLSF (биты19-26). Все эти прогнозные leap-second поля теряются, блок WN/WNLSF искажается. Все сиблинг-кейсы (MT30/31/37) маппят слова в последовательные слоты — MT33 единственный с дублем и нулевым data[8].

**Фикс:** `data[7] = self.utc_message[3];` (и `data[8] = 0`).

---

### H13. B-CNAV3 `gf6_int_mul`: identity-антилог и испорченный лог-стол ломают LDPC-кодирование
**Файл:** `src/bcnav3bit.rs:290-302`

`gf6_int_mul` вычисляет `E2V_TABLE[(V2E_TABLE[a]+V2E_TABLE[b]) % 63]`, но `E2V_TABLE` — это identity `[0,1,2,...,63]`, а `V2E_TABLE` не биекция (`V2E[9]==V2E[14]==15`, `V2E[18]==V2E[28]==16`; всего 21 дубль). Численно подтверждено: таблицы не имеют мультипликативной единицы, т.е. не задают поле; 3882 из 3969 ненулевых произведений отличаются от корректных. Корректные таблицы GF(2⁶) уже существуют в `bcnavbit.rs:132-146`. `get_frame_data → ldpc_encode (line 344) → gf6_int_mul (line 266)` — живой путь, результат пишется в nav_bits (352-354), так что каждый символ паритета B2b (B-CNAV3) LDPC испорчен.

**Фикс:** заменить таблицы на проверенные `BCNavBit::E2V_TABLE/V2E_TABLE` (с учётом конвенции индексации того файла) или вызывать общий `BCNavBit::gf6_int_mul`.

---

### H14. B-CNAV2 `append_word`: читает MSB (бит 32) и неверно продвигает индекс источника
**Файл:** `src/bcnav2bit.rs:338-352`

Два дефекта: (1) строка 338 извлекает `(source[src_index] >> (32 - bits_to_copy)) & mask`, читая старшие биты 32-битного слова, но payload упакован в низкие 24 бита (документировано `bcnav2bit.rs:183`; эталон `bcnavbit.rs:415` сдвигает `>> (remain_bits - fill_bits)` от 24). Даже полный 24-битный чанк даёт сдвиг `32-24=8` → `(source>>8)`, теряя низкие 8 бит payload. (2) `src_index` продвигается только при `bits_to_copy >= 24`, поэтому при разбиении на под-24-битные чанки слово рассинхронизируется по выравниванию. Трассировка реальных вызовов (msg 10/11/30/31/32/33/34/40) показывает, что цикл завершается и хвостовые слова НЕ теряются (формулировка исходной заявки «бесконечно перечитывается» неточна), но payload искажён по битам. B-CNAV2 (B2a) — живой формат (`ifdatagen.rs:586,1145`); эфемерис1/2, clock, almanac получают неверные данные.

**Фикс:** портировать логику `bcnavbit::append_word` (24-битное окно источника, счётчик remaining-bits на слово, продвигающий src_index при исчерпании), либо делегировать в общую реализацию.

---

### H15. GLONASS-эфемерида: знаковый бит X/Y/Z потерян (кодируется только магнитуда)
**Файл:** `src/gnavbit.rs:246-248, 263-265, 283-285`

`ComposeStringEph` кодирует `uint_value = (eph.x.abs() / 2^-11).round() as u32`, упаковывая только магнитуду в 27-битные поля, и НИКОГДА не ставит знаковый бит. ICD GLONASS (Ed. 5.1 Table 4.5) задаёт координаты как sign-magnitude (1 знак + 26 магнитуда). Almanac-путь в этом же файле (`ComposeStringAlm:333-369`) корректно ставит знаки (`|= if clock_error < 0.0 { 1 << 9 }`), доказывая намерение. GLONASS ECEF-координаты отрицательны примерно на половине орбиты по каждой оси — приёмник декодирует |x| и помещает спутник в неверный октант.

**Фикс:** ставить MSB (бит 26) 27-битных полей при отрицательной координате, как для almanac.

---

### H16. GLONASS-эфемерида: скорость/ускорение/часы кодируются two's-complement вместо sign-magnitude
**Файл:** `src/gnavbit.rs:250-255, 267-272, 287-292, 297-308`

vx/vy/vz, ax/ay/az, tn/dtn/gamma кодируются `int_value = (value/scale).round() as i32` и `COMPOSE_BITS!(int_value,...)`, где `COMPOSE_BITS` делает `($value as u32) & mask` — для отрицательного i32 сохраняется two's-complement битпаттерн. ICD GLONASS задаёт sign-magnitude (знак в MSB поля + абсолютная магнитуда). Для 5-битного поля accel с int_value=-1: two's-complement = 0b11111 (31), sign-magnitude = 0b10001 (17). Эти величины физически регулярно отрицательны (компоненты скорости, clock bias tn). Almanac-путь делает sign-magnitude корректно, подтверждая нужную конвенцию. Тест `nav_message_formation.rs:277` зануляет vx/vy/vz/ax/ay/az, так что баг не ловится.

**Фикс:** кодировать как магнитуда + явный знаковый бит в MSB поля, зеркаля `ComposeStringAlm`.

---

### H17. GLONASS tk: кодируется плоский tk/30 вместо битового layout часы|минуты|30с
**Файл:** `src/gnavbit.rs:244`

Строка 244: `COMPOSE_BITS!(ephemeris.tk / 30, 2, 12)`. **Важная поправка к исходной заявке:** парсеры RINEX (`json_interpreter.rs:2763, 3402`) хранят tk УЖЕ в битовом layout: `tk = ((hours<<7)|(minutes<<1)|half_minutes)`. Поэтому деление этого упакованного значения на 30 портит его целиком: для tk_seconds=3630 хранится tk=129, кодер пишет 129/30=4, а должно остаться 129. Сиблинг tb хранится в сырых секундах (`json_interpreter.rs:3405`) и кодер корректно делит на 900 (`gnavbit.rs:261`); tk же надо писать напрямую. tk искажён практически для всех реальных GLONASS-эфемерид. Тесты используют синтетический tk:1200 и не декодируют поле.

**Фикс (корректный для текущей конвенции хранения):** `COMPOSE_BITS!(ephemeris.tk, 2, 12)` — tk уже несёт layout `(h<<7)|(m<<1)|s30`. (Предложенный в заявке пересчёт из секунд привёл бы к двойной упаковке.)

---

### H18. Galileo I-NAV word 5: ионосфера NeQuick (ai0/ai1/ai2) молча затирается в `compose_eph_words`
**Файл:** `src/inavbit.rs:768-769`

`ephdata[16] &= 0x03ffffff;` (сохранить низкие 26 бит — поле ионосферы, записанное `set_iono_utc`) немедленно перезаписывается `ephdata[16] = 0x14000000;` (плейн-присваивание, не OR), отбрасывая маскированное значение — классический dead store. Маска `0x03ffffff` иначе бессмысленна. Кроме того, строка 773 `ephdata[17] = COMPOSE_BITS!(...)` так же затирает `ionowords[1]`. В сигнальном пути `set_iono_utc` (`ifdatagen.rs:1518`) выполняется до per-SV `set_ephemeris` (`1610`), поэтому передаваемые ai0/ai1/ai2 в I-NAV word 5 всегда нули. Тесты I-NAV не вызывают `set_iono_utc`, баг не ловится.

**Фикс:** `ephdata[16] = (ephdata[16] & 0x03ffffff) | 0x14000000;` (сохранить ионосферу, выставив type=5 в верхних 6 битах). Аналогично для строки 773. Добавить тест: загрузить iono-параметры, построить word 5, декодировать ai0/ai1/ai2 != 0.

---

### H19. Горизонтальная круговая траектория: вращение скорости и смещение позиции в несовместимых системах отсчёта
**Файл:** `src/trajectory.rs:673-701`

В `TrajectoryHorizontalCircular::get_pos_vel` скорость вращается вокруг сырой ECEF-оси z (`new_vx = vx·cos - vy·sin`), а смещение позиции вычисляется в heading-выровненной системе (`dx = radius·sinθ`, `dy = radius·(1-cosθ)`) и добавляется прямо к ECEF x/y. (1) При t=0 производная позиции = `(speed,0,0)`, а выходная скорость = реальная `(vx,vy,vz)` — расходятся, если начальная скорость не строго вдоль +ECEF-x; интегрированный трек и сообщаемая скорость взаимно противоречивы. (2) Вращение ECEF-скорости вокруг ECEF-z — горизонтальный поворот только на полюсе; на любой широте ECEF-z ≠ локальная вертикаль. Сиблинг `TrajectoryVerticalAcc:413-418` корректно выводит локальную вертикаль через `calc_up_vector(ecef_to_lla(...))`, доказывая, что фреймворк знает про различие. Смещение `(dx,dy)` никогда не вращается в начальный курс и не проходит через имеющийся ENU `convert_matrix`. Путь не покрыт тестами.

**Фикс:** строить круговое движение в ENU через `base.convert_matrix` и `base.get_speed_projection` (начальный курс), вращать смещение `(radius·sinθ, radius·(1-cosθ))` в направление начального курса, затем отображать в ECEF — как делают остальные типы сегментов.

---

### H20. Test gap: единственный Galileo FEC-декод тест `#[ignore]`'нут; F/NAV FEC и CRC вообще не тестируются
**Файл:** `tests/signal_quality.rs:1034-1035`

`test_galileo_inav_viterbi_decode` — единственный тест, который деинтерливит, Viterbi-декодит и CRC24Q-проверяет реальный I-NAV битстрим — помечен `#[ignore]` и known-failing (0/15 страниц CRC OK, even-part data mismatch на всех страницах). F/NAV не имеет FEC/CRC-декод теста вовсе (только проверка 12-битного sync и плотности битов). Из-за этого CRC24Q-несовпадение, F/NAV G1-ошибка (C1) и отсутствие инверсии G2 проходят CI незамеченными. Round-trip декодеры в `nav_decode.rs` работают на сырых word-массивах и структурно не касаются свёрточного FEC.

**Фикс:** исправить CRC24Q и even-part reconstruction, снять `#[ignore]` с I-NAV Viterbi-теста; добавить аналогичный F/NAV тест deinterleave+Viterbi+CRC с корректными полиномами (G1, G2 инвертирован).

---

## Системные темы

Повторяющиеся первопричины, проходящие через находки:

1. **Ошибки модели фазы / модуло на границах мс (H1, H2, H3).** Несущая теряет дробный IF-вклад каждую мс; code phase обёртывается по неверной длине (`pilot_length=1` для L2P, `0.0` center freq для G3). Общий паттерн — кусочно-помиллисекундная аккумуляция с реанкорингом, теряющая величины, которые эталонный непрерывный генератор `gps_pilot/mod.rs` сохраняет. Несущая и код должны аккумулироваться непрерывно с обёрткой по истинному периоду сигнала.

2. **Дивергентные дубль-реализации одного алгоритма (H4, H5–H11, H13, H14).** Существует «правильная» эталонная версия (`coordinate::glonass_sat_pos_speed_eph`, `cnavbit.rs`, `bcnavbit.rs`, `inavbit::convolution_encode`), а параллельная копия (`satellite_param.rs`, `l5cnavbit.rs`, `bcnav3bit.rs`, `bcnav2bit.rs`, `fnavbit.rs`) отклоняется в масштабах, полиномах, GF-таблицах или конвенции единиц. Дублирование кода — главный источник High-багов.

3. **Неверная конвенция единиц / масштаба vs ICD (H4, H6, H9, H10, H11, H17).** Лишний ×1000 (км/м), дважды учтённые секунды, степени двойки CNAV (2⁻³⁴ vs 2⁻³⁵), sqrtA vs ΔA, плоский tk vs битовый layout — все из-за рассинхрона между парсером и потребителем по единицам/конвенции хранения.

4. **Потеря/отсутствие инициализации полей (H8, H18).** GPS-контейнер используется для BeiDou без инициализации `tgd2`; ионосфера I-NAV word 5 затирается dead store. Lossy-каст `f64→u16 as` молча портит iodc.

5. **Sign-magnitude vs two's-complement в кодировании nav-сообщений (H12, H15, H16).** Знаковые поля GLONASS/CNAV кодируются неверным представлением; almanac-путь делает правильно, эфемерис-путь — нет (несогласованность внутри одного файла).

6. **Несогласованность видимость↔генерация (H5).** Набор видимых спутников и набор генерируемых считаются разными позиционными моделями (заглушка vs RK4).

7. **Системный test gap для nav-данных (H20 и сопутствующие).** Существующие тесты проверяют только sync-паттерны, CRC-наличие, плотность битов и round-trip против тех же нестандартных хелперов; единственный реалистичный full-pipeline FEC/CRC тест отключён. Поэтому весь класс ICD-багов (C1, H9–H18) проходит CI.

---

## Рекомендации

**Чинить в первую очередь (ломает сигнал/время для shipped-пресетов):**
1. **C1** — F/NAV G1/G2 (заменить ОБЕ маски на `0x4F` и `0x6D^1`); E5a иначе недекодируем.
2. **H1** — аккумуляция полной фазы (Doppler+IF) в якорь; ломает PLL для `quad_l1g1.json` и всего GLONASS FDMA.
3. **H2, H3** — G3 center freq и L2P wrap; `glo_g3.json`/`gps_l2p.json` выдают полностью недетектируемый сигнал.
4. **H4, H5** — GLONASS-пропагатор ×1000 и no-op заглушка видимости; делегировать обе в `coordinate::glonass_sat_pos_speed_eph` (убивает дубликаты).
5. **H7** — RINEX `>` → `>=`; триггерится репозиторной фикстурой, вероятно объясняет недосчёт Galileo.

**Во вторую очередь (ICD-некорректные nav-данные, ломают реальные приёмники):**
6. **H8** — BeiDou TGD2/IODC.
7. **H9–H12** — переписать L5 CNAV MT10/MT11/MT30, опираясь на `cnavbit.rs`; и MT33 UTC.
8. **H13, H14** — B-CNAV3 GF-таблицы и B-CNAV2 `append_word` (переиспользовать `bcnavbit`).
9. **H15–H17** — GLONASS sign-magnitude для X/Y/Z и v/a/clock, фикс tk.
10. **H18** — dead store ионосферы I-NAV word 5.

**В третью очередь:**
11. **H6** — `utc_to_gps_time` (латентно при целых секундах, но корруптит всю базу времени).
12. **H19** — горизонтальная круговая траектория через ENU.

**Тесты, которые нужно добавить (закрывают системный gap #7):**
- Снять `#[ignore]` с `test_galileo_inav_viterbi_decode` после фикса CRC24Q; добавить F/NAV deinterleave+Viterbi+CRC тест с G1=0x4F/G2 инверсия (ловит C1, H20).
- Bit-true decode-тесты для L5 CNAV MT10/MT11/MT30 (масштабы af0/af1/af2, ΔA, TGD/ISC) — ловят H9–H11.
- Round-trip BeiDou D1/D2: проверять tgd2 != 0 после парсинга реального RINEX (H8).
- Тест `utc_to_gps_time` с дробной стартовой секундой (напр. 30.5), assert `SubMilliSeconds < 1.0` и round-trip (H6).
- RINEX-парсинг тест на right-trimmed строках длиной 23/42/61, assert `eph.week != 0` (H7).
- GLONASS nav decode: декодировать X/Y/Z/v/a/tk из эмитированных бит и сравнить с эфемеридой, включая отрицательные координаты (H15–H17).
- I-NAV word 5: загрузить iono-параметры, построить word 5, декодировать ai0/ai1/ai2 != 0 (H18).
- GLONASS визибилити-vs-генерация: сравнить позиции из обоих путей при |t−tb| ~ 600с (H5).

Ключевые затронутые файлы: `src/sat_if_signal.rs`, `src/fnavbit.rs`, `src/satellite_param.rs`, `src/ifdatagen.rs`, `src/gnsstime.rs`, `src/json_interpreter.rs`, `src/cnavbit.rs`, `src/l5cnavbit.rs`, `src/bcnav3bit.rs`, `src/bcnav2bit.rs`, `src/gnavbit.rs`, `src/inavbit.rs`, `src/trajectory.rs`, `tests/signal_quality.rs`.


---

## Medium (детально)

- **src/satellite_param.rs:484** — BeiDou B2AB Doppler uses wrong (B1C) wavelength
  - get_doppler() (line 570-572) divides RelativeSpeed by the wavelength to produce the Doppler that drives carrier phase. For B2AB the carrier is generated at ~1191.795 MHz but Doppler is scaled by the B1C wavelength (1575.42 MHz) — about a 1.32x error in Doppler magnitude. A receiver tuned to B2AB would see Doppler ~32% too large, breaking acquisition/tracking. The center-frequency table and wavelen
  - *Фикс:* Add `SIGNAL_INDEX_B2AB => LIGHT_SPEED / FREQ_BDS_B2AB,` to the BdsSystem match in get_wave_length(), consistent with the center-frequency tables.
- **src/ifdatagen.rs:1206, 1400, 1440 (and 4099)** — Ephemeris epoch selection uses UTC-scale time (leap=false) while signal generation uses GPS-scale time (leap=true) — 18 s mismatch
  - src/gnsstime.rs:209-222 shows use_leap_second=true adds the leap seconds to reach GPS time; false yields a UTC-scale value 18 s earlier. tests/time_system_baseline.rs:22-24 fixes the canonical GPS time for 2025-06-05 at MilliSeconds=381_948_250 (leap applied). select_per_satellite_and_fill (src/ifdatagen.rs:1400 GPS, 1440 GAL) and select_global_epochs_and_fill (1206) use leap=false → 381_930_250, 
  - *Фикс:* Use utc_to_gps_time(utc_time, true) consistently for GPS and Galileo epoch-selection targets (lines 1206, 1400, 1440) and in src/json_interpreter.rs:1149 and is_ephemeris_within_time_window (utc_to_gps_time(*target_time, false)), so the selection/filter basis matches self.cur_time (leap=true) used i
- **src/ifdatagen.rs:1459-1477** — Active GLONASS per-satellite epoch selection ignores Day, matching only seconds-of-day within ±1800 s
  - GLONASS tb is a time-of-day reference (0..86400 s); the same tb recurs every day. The companion routine find_glo_ephemeris (src/ifdatagen.rs:243-253) explicitly adds day_diff*86400 to disambiguate days, proving the Day component is needed for correct matching, but the active selection path here drops it. The comment acknowledges 'many RINEX do not contain a correct day for GLONASS', yet for files 
  - *Фикс:* Incorporate eph.day (cumulative day from the GLONASS epoch in the corrected GLONASS time) into the diff, e.g. diff = (gtime.Day - eph.day)*86400 + (req_sec - eph.tb), mirroring find_glo_ephemeris, falling back to seconds-of-day only when the RINEX lacks day info.
- **tests/time_system_baseline.rs:19-57** — No test verifies ephemeris-selection time basis matches the generation clock
  - The leap=false vs leap=true split between selection (src/ifdatagen.rs:1400/1440, src/json_interpreter.rs:1149) and generation (src/ifdatagen.rs:4099) is exactly the class of bug that produced the historical '4-year offset / satellites below horizon' regression documented in CLAUDE.md. The existing tests pass with the flags hard-coded per call, so they cannot catch this cross-module inconsistency.
  - *Фикс:* Add a test that, for a fixed preset UTC, asserts the epoch-selection target seconds-of-week equals self.cur_time's seconds-of-week (i.e. both leap-true), or directly asserts the delta_t a satellite sees during generation equals the |Δt| reported by epoch selection within <1 s.
- **src/satellite_signal.rs:437** — Pilot secondary-code phase uses unadjusted transmit_time while data channel uses leap-second-adjusted time (BeiDou misalignment)
  - For BeiDou B1C (code_length=10, secondary_code_length=1800) and B2a (code_length=20, secondary_code_length=100), the data bits are taken at transmit_time-14000 ms while the pilot secondary code is taken at transmit_time, a 14000 ms (=1400 secondary positions for B1C, 700 for B2a) offset between the two channels of the same satellite. The data and pilot of one SV are transmitted on the same code ep
  - *Фикс:* Use the same time base as the data channel for the secondary-code index: replace `transmit_time.MilliSeconds` at line 437 with `transmit_time_adj.MilliSeconds` (matching how frame_number/bit_number are derived).
- **src/sat_if_signal.rs:740-744** — QMBOC/TMBOC BOC(6,1) chip selection driven by millisecond index instead of spreading-chip index
  - TMBOC/QMBOC (per IS-GPS-800 / BDS-SIS-ICD) places the BOC(6,1) component on specific spreading-code CHIP positions (4 out of every 33 spreading symbols), not on millisecond indices. The 33 distinct symbol_pos values here cycle once per 10 ms and are constant within a code period, so the intra-code-period BOC(6,1)/BOC(1,1) chip pattern required by the ICD is never produced; the BOC(6,1) subcarrier 
  - *Фикс:* Derive the symbol index from the spreading chip number within the code period (e.g. from pilot_chip / and its position modulo 33) rather than from MilliSeconds, and apply the BOC(6,1) sub-chip flip only on the ICD-specified chip positions; the BOC(6,1) sub-chip square wave (chip_raw.rem_euclid(12) >
- **tests/signal_assembly_baseline.rs:202-232** — No test covers carrier-phase continuity for sat_if_signal with a fractional (non-kHz) IF frequency
  - The critical carrier-phase discontinuity (this report's finding #1) only manifests when if_freq is not a multiple of 1000 Hz and only across consecutive milliseconds. A single-ms energy check is structurally incapable of catching it, which is why the shipped quad_l1g1.json preset (0.5 cycle/ms jump on GPS/BDS/GAL) is not flagged by the suite. This is a high-value gap given the file is the project'
  - *Фикс:* Add a test that runs get_if_sample_cached for many consecutive milliseconds with a deliberately fractional if_freq (e.g. 562500 Hz, frac=0.5), reconstructs the per-sample phase, and asserts the phase delta between the last sample of ms N and the first sample of ms N+1 is within one phase_step (max j
- **src/satellite_param.rs:344 (and 256-385)** — Relativistic clock correction is a permanent no-op (eph.Ek never populated) in get_satellite_param
  - The relativistic eccentricity correction Δt_r = F·e·√A·sin(E) (F=-4.442807633e-10 s/√m, GPS ICD-200) is mandatory and reaches ±20-45 ns (±6-13 m pseudorange) at e~0.02. ifdatagen.rs contains NO independent relativistic term (grep for 4.442807633/WGS_F_GTR/Ek in ifdatagen returns nothing), so TravelTime from get_satellite_param is the only place it could be applied — and it evaluates to zero becaus
  - *Фикс:* Make the local gps_sat_pos_speed_eph write Ek/Ek_dot back (take `&mut GpsEphemeris` or return ek), or compute the eccentric anomaly locally inside get_satellite_param and use that value in the relativistic term instead of `eph.Ek`. Reuse `crate::coordinate::gps_sat_pos_speed_eph` which already popul
- **src/satellite_param.rs:344 vs 1716-1719** — Relativistic correction sign inconsistent between get_satellite_param and SatelliteParamCalculator
  - The GPS/ICD relativistic correction is Δt_r = F·e·√A·sin(E) with F = -4.442807633e-10 (negative). The clock correction is subtracted from travel_time, and dt_r adds to clock error, so the term should be `travel_time -= WGS_F_GTR*...` (as in line 1716). Line 344 uses the positive literal, i.e. the wrong sign; only the always-zero Ek (finding above) currently masks the resulting ~±2×13 m error. If E
  - *Фикс:* Replace the positive literal at line 344 with `WGS_F_GTR` (the negative constant) so it matches line 1716 and the ICD sign convention.
- **src/coordinate.rs:1001-1016** — cis_to_cts writes rotated POSITION into the acceleration output array instead of acceleration
  - KinematicInfo has no acceleration members, so there is no acceleration to rotate; the code substitutes position. Any caller requesting GLONASS acceleration receives the rotated position vector (~2.6e7 magnitude) instead of an acceleration (~0.6 m/s^2), an error of ~7 orders of magnitude. Currently latent/dead because every caller of glonass_sat_pos_speed_eph in the focus files passes acc=None (coo
  - *Фикс:* Either remove the acceleration branch from cis_to_cts (and the acc plumbing through glonass_sat_pos_speed_eph) since KinematicInfo cannot carry acceleration, or extend the state to include acceleration components (ax,ay,az propagated by RK4) and rotate the true acceleration. At minimum, do not emit 
- **tests/rinex_parser_baseline.rs:75-117, 43-73** — RINEX parser tests verify counts/SVIDs but never orbital/TGD field values
  - Because no test pins the actual BeiDou TGD2 / IODC values or GPS/Galileo orbital fields to the fixture, the BeiDou TGD2=0 bug (data[26] miscast to iodc) and the field off-by-one boundary issues go undetected. A regression that zeroes or swaps an ephemeris field would still pass the current suite.
  - *Фикс:* Add assertions that parse a known fixture record and compare decoded fields to hand-verified values: e.g. for C01 in rinex_v3_20251560000.rnx assert tgd1≈-5.2e-9, tgd2≈-9.7e-9, sqrtA≈6493.478..., toe=345600; for G01 assert sqrtA≈5153.7217, ecc≈5.6157e-4, week=2369.
- **src/l5cnavbit.rs:81, 219, 508-509** — L5 CNAV almanac toa stored in u8, truncating the 13-bit WNa and 8-bit toa to 8 bits total
  - The CNAV almanac time-of-applicability combines WNa (13 bits) and toa (8 bits) = 21 bits, as the L2C path correctly stores in a u32 (cnavbit.rs:44,219). Storing it in a u8 truncates `(week<<8)+(toa>>12)` to the low 8 bits, discarding WNa entirely. Then get_message_payload does `self.toa as u32 >> 13`, but a u8 value (<=255) shifted right by 13 is always 0, so WNa is always transmitted as zero, and
  - *Фикс:* Change the field to `toa: u32` (as in cnavbit.rs), store the full 21-bit (WNa<<8)|toa value, and split it correctly in get_message_payload.
- **src/l5cnavbit.rs:151-163** — L5 CNAV CRC excludes preamble+PRN (computed over 262 of 276 bits), diverging from CNAV/ICD
  - Per IS-GPS-705, the L5 CNAV CRC-24Q is computed over the full 276-bit message (preamble + PRN + message type + TOW + alert + data). The in-repo L2C reference cnavbit.rs includes the preamble/SVID header in encode_data[0] and CRCs all 276 bits (cnavbit.rs:289,137). The L5 path seeds the CRC at the message-type field, so the transmitted CRC does not protect the preamble/PRN and will not match a conf
  - *Фикс:* Include the 8-bit preamble (0x8B) and 6-bit PRN in the CRC input block and compute CRC over the full 276 bits, as cnavbit.rs does.
- **tests/icd_conformance_fixtures.rs:283-324** — No round-trip / scale-factor test coverage for L2C CNAV, L5 CNAV, and CNAV-2 ephemeris encoding
  - The LNAV path has decode_lnav_stream123 + compare_with_original round-trip tests that would catch a wrong scale factor or bit offset. CNAV, L5 CNAV and CNAV-2 have no equivalent decoder/round-trip check, which is precisely why the L5 scale-factor errors (af0/af1/af2, M0/omega/i0/w, TGD) and the cnavbit MT33 data[7] duplication go undetected. The L5 golden-snapshot assertions (l5_bits[..48], [552..
  - *Фикс:* Add decode-and-compare round-trip tests for CNAV MT10/11/30 and L5 CNAV MT10/11/30 against the source ephemeris (LSB/2 tolerances), as done for LNAV in src/lnavbit.rs tests via nav_decode.
- **src/bcnavbit.rs:269** — TGD_B2ap field encoded from wrong source array element (tgd_ext[1] instead of tgd_ext[2])
  - Line 266-270 establishes the intended mapping: isc_b1c=tgd_ext[0]-tgd_ext[1], tgd_b1c=tgd_ext[1], isc_b2a=tgd_ext[2]-tgd_ext[3], tgd_b2a=tgd_ext[1]. The B2a TGD must come from the B2a slot (tgd_ext[2]), not the B1C slot (tgd_ext[1]). Per BDS-SIS-ICD-B2a, TGD_B2ap is a distinct group-delay parameter. Using tgd_ext[1] makes the transmitted TGD_B2ap equal to TGD_B1Cp, which is physically wrong and wi
  - *Фикс:* Change line 269 to `let tgd_b2a_scaled = Self::unscale_int(eph.tgd_ext[2], -34);`
- **src/inavbit.rs:667-672, 768-773** — set_iono_utc ionosphere data in word 5 is wiped if set_ephemeris runs afterward
  - If the caller invokes set_iono_utc before set_ephemeris (a natural order), the iono parameters placed into words [16]/[17] are destroyed by the unconditional assignments in compose_eph_words. The two writers both own word-5 bits but one clobbers the other, producing pages with zero/garbage ionosphere data depending on call order.
  - *Фикс:* Make compose_eph_words use |= for the word-5 region (preserving any pre-set iono bits) or re-apply iono after composing ephemeris, and document the required call order.
- **src/prngenerate.rs:408-411** — GLONASS G3 (L3OC) Gold code is degenerate: G1 register outputs a constant, code repeats with period 341
  - Simulated with exact Rust LsfrSequence semantics (u32 state, feedback masked by polynomial, output = bit depth-1): G1 first 1023 outputs are all 1 (sum=1023). The resulting code has measured period 341, not 10230. The prn_generation_baseline test passes only because balance_ratio (0.062) is below its 0.10 threshold and it never checks periodicity/cross-correlation. The CLAUDE.md claims 'G3/L3OC us
  - *Фикс:* Use the correct GLONASS L3OC primary-code generator (the real L3OC code is generated from longer registers / a Kasami-like construction, not two 10-bit registers). At minimum the G1 polynomial/init must produce a genuine maximal-length sequence and the register width must support a 10230-chip code. 
- **src/avx512_intrinsics.rs:160-167** — Write through *const-derived pointer into immutable array is UB in fast_sin_cos_avx512
  - Mutating data behind a shared/`*const` reference to an immutable binding without `UnsafeCell` is undefined behavior in Rust. The optimizer is entitled to assume `indices_array` never changes and may replace the reads at lines 164-167 with the initializer value 0, so every lookup returns `sin_lut[0]`/`cos_lut[0]` (= 0 and 1) instead of the gathered values. Even where it happens to work today, it is
  - *Фикс:* Declare `let mut indices_array = [0i32; 16];` and store via `_mm512_storeu_si512(indices_array.as_mut_ptr() as *mut __m512i, masked_indices);`.
- **src/avx512_intrinsics.rs:98-128, 238-250, 254-270** — AVX-512 functions lack #[cfg(target_arch = "x86_64")] gating and fail to compile on non-x86 targets
  - The `core::arch::x86_64::*` import (line 5) is `#[cfg(target_arch = "x86_64")]`, so on any non-x86_64 target the `_mm512_*` calls in these three bodies are unresolved, and `enable = "avx512f"` is not a valid target feature for those architectures. `cargo build --target aarch64-...` (Android/ARM targets are installed, and CLAUDE.md references phone-receiver testing) therefore fails to compile. The 
  - *Фикс:* Add `#[cfg(target_arch = "x86_64")]` above the `#[target_feature(enable = "avx512f")]` on all three functions, matching the rest of the module.
- **src/almanac.rs:40-70** — check_almanac_type can never return AlmanacBds, so BeiDou almanac files are misclassified and read_almanac_bds is unreachable
  - read_almanac_bds (line 254) and get_almanac_bds (line 271) parse a whitespace-separated BeiDou almanac with >=12 numeric columns and no leading '*' or '<'. A real BeiDou almanac file fed to read_almanac will instead fall into the else-branch of check_almanac_type: parts[1] is the second numeric column, which is not in DD.MM.YYYY form, so it returns AlmanacUnknown and read_almanac returns 0, silent
  - *Фикс:* Add BeiDou detection to check_almanac_type (e.g. detect a numeric first token in the valid SVID range with the expected column count, or a BDS-specific header/marker), returning AlmanacType::AlmanacBds so the existing read_almanac_bds path is reachable.
- **src/galileo_pilot/mod.rs:85-93** — Galileo E1-C CS25 secondary code extracted MSB-first — reversed/wrong sequence vs. rest of codebase
  - The documented CS25 bit string `0011100000001010110110010` is reproduced exactly by LSB-first extraction `(0x9b501c >> i) & 1` (verified: this yields the documented string), which is what satellite_signal.rs and gnss_pilot.rs both use. The galileo_pilot MSB-first formula `(0x9b501c >> (24-i)) & 1` instead yields `0100110110101000000011100` — the time-reversed sequence, a completely different 25-ch
  - *Фикс:* Change the extraction to LSB-first to match satellite_signal.rs and gnss_pilot.rs: `if (SECONDARY_CODE_RAW >> idx) & 1 != 0 { -1.0 } else { 1.0 }` (drop the `24 -`). Remove the unused/misleading SECONDARY_CODE = 0x0E64C2E0 constant and fix the `MSB first` comments. Add a unit test asserting the firs
- **src/galileo_pilot/mod.rs:85-93** — Galileo pilot secondary-code path has no test coverage
  - The MSB-vs-LSB secondary-code defect (separate finding) shipped precisely because nothing asserts the emitted CS25 sequence. A direct test comparing the first 25 secondary bits against either the documented string `0011100000001010110110010` or against gnss_pilot.rs's `(CS25 >> idx) & 1` extraction would have caught the reversal immediately. Lack of any acquisition test for galileo_pilot also mean
  - *Фикс:* Add a unit test in src/galileo_pilot/mod.rs asserting get_secondary_code_bit(0..25) matches the documented CS25 sequence and matches gnss_pilot's extraction, plus an acquisition-style integration test (mirroring tests/gps_pilot_test.rs) that generates a short Galileo E1 signal and confirms the expec
- **src/coordinate.rs:463-486** — GLONASS RK4 propagator drops the residual (sub-30s) integration step on the tb path
  - The integer-truncated step count leaves up to 30 s of orbital motion un-propagated (up to ~108 km at 3.6 km/s) on the first (flag&2==0) call, while the CIS→CTS rotation uses the full delta_t, so the position and the frame-rotation correspond to different times. The alternate tc-prediction branch (line 508) DOES apply the residual via a single runge_kutta(delta_t1), proving the tb branch is inconsi
  - *Фикс:* After the coarse loop, apply the residual step, e.g. `let residual = delta_t - (step_number as f64)*COARSE_STEP; if residual != 0.0 { runge_kutta(residual, &mut state); }`, mirroring the tc-prediction branch, before writing eph.PosVelT and calling cis_to_cts.
- **src/ifdatagen.rs:762, 2234, 2531** — IQ4 quantization scale (3.0) leaves the AGC-targeted signal at ~1 of 7 magnitude levels
  - For the IQ4 output path the 4-bit dynamic range is grossly under-utilized: most samples round to magnitude 0 or 1, effectively ~1-bit quantization, destroying SNR for IQ4 captures. The IQ8 path is fine because QUANT_SCALE_IQ8 = 127 maps RMS 0.25 to ~32 counts. The single shared target_rms=0.25 is calibrated for the 127-scale path, not the 3.0-scale IQ4 path; either the scale or the target RMS must
  - *Фикс:* Pick QUANT_SCALE_IQ4 so that ~3σ reaches near full scale, e.g. scale ≈ 7/(3*target_rms) ≈ 9.3 (or use a format-specific AGC target). Verify against quantize_samples_iq4 test expectations (the existing test uses inputs |value|>=1 so it may need updating).
- **src/sat_if_signal.rs:496-501 and 656-669** — Carrier phase accumulator omits IF-frequency term -> per-millisecond phase discontinuity for off-center signals
  - Working modulo 1 cycle, the end-of-ms phase is `-scp_k + doppler/1000 + if_freq/1000` while the start-of-next-ms phase is `-scp_{k+1} = -scp_k + doppler/1000`; the difference is `frac(if_freq/1000)` cycles, a discontinuity at every ms boundary whenever if_freq is not a whole multiple of 1000 Hz. if_freq is set as `(signal_center_freq - output_center_freq) as i32` (sat_if_signal.rs:316/344); it is 
  - *Фикс:* Accumulate the full carrier phase (Doppler + IF) into start_carrier_phase across milliseconds: e.g. `self.start_carrier_phase -= doppler_cycles_per_ms + self.if_freq as f64 / 1000.0;` (matching the per-sample phase_step), or maintain a separate continuous IF phase. Add a carrier-phase-continuity tes
- **tests/signal_quality.rs:373-435** — test_galileo_inav_all_svs_crc verifies only the sync pattern, never CRC or page content
  - The sync pattern is a hardcoded literal copied straight from inavbit.rs:50 (SYNC_PATTERN) into nav_bits at get_frame_data lines 583-586/594-597; asserting it round-trips proves nothing about the encoder. Running the ignored Viterbi test now shows CRC MISMATCH and 'even=false' data mismatch for every page tested (tow=1..29), so the entire Galileo I-NAV page content is effectively unverified by the 
  - *Фикс:* Either rename this test to reflect that it only checks sync, or fold in the de-interleave + Viterbi + CRC24Q logic from the ignored test and assert CRC passes. Investigate and fix the even-half data mismatch surfaced by test_galileo_inav_viterbi_decode (tail bits pass and the odd half matches, isola
- **tests/icd_conformance_fixtures.rs:86-94, 418-449** — ICD conformance golden values for BeiDou B1C are captured from buggy/non-representative encoder output
  - flag=255 is not a valid B-CNAV1 SatType (ICD: 1=GEO,2=IGSO,3=MEO; 0 reserved), so the fixture exercises a code path that cannot occur with correctly converted data (to_gps_ephemeris sets flag=sat_type+1). The golden LDPC vectors therefore lock in encoder output for an invalid input and validate self-consistency, not ICD conformance; if the axis-reference logic is later fixed, these 'bit-exact' ass
  - *Фикс:* Drive the fixture through BeiDouEphemeris::to_gps_ephemeris() (or set flag to a valid 1/2/3 consistent with the chosen axis/SVID) and recompute golden vectors from a known-correct reference (e.g. cross-checked against the BDS-SIS-ICD-B1C test vectors), not from current output. Add separate cases for
- **src/bcnav1bit.rs:658-669** — BeiDou B-CNAV1 SatType field encoded as 0/1/2 (GEO=0) instead of ICD 1/2/3, inconsistent with axis-reference logic
  - Per BDS-SIS-ICD-B1C the SatType code is 01=GEO, 10=IGSO, 11=MEO; encoding GEO as 0 (a reserved value) is non-conformant, and a real B1C receiver decoding SatType=0 would treat it as reserved/invalid. The mismatch with bcnavbit.rs's 1/2/3 convention also means the field broadcast in Frame 2 and the axis-reference used in Ephemeris1 describe the satellite type with two different numbering schemes wi
  - *Фикс:* Encode SatType per ICD as 1 (GEO), 2 (IGSO), 3 (MEO), consistent with eph.flag used for the reference-axis selection; alternatively derive both from a single source of truth (eph.flag).
- **tests/upstream_fixes.rs:200-208** — Placeholder test test_find_ephemeris_equal_time_diff asserts nothing (assert!(true))
  - The behavior it claims to protect is real and load-bearing: IFDataGen::find_ephemeris (ifdatagen.rs:185/201/217) uses `abs_diff <= best_time_diff`, so a later equal-diff ephemeris replaces an earlier one. A regression to `<` would silently change which ephemeris every satellite uses, and this test would still pass. (Note also the divergent duplicate in json_interpreter.rs:1170/1191/etc., which use
  - *Фикс:* Construct two GpsEphemeris with equal |toe - target| at different iteration positions, call IFDataGen::find_ephemeris, and assert the later one's distinguishing field is returned. Reconcile the json_interpreter find-ephemeris tie-break (< vs <=) with IFDataGen so both pick the same epoch.
- **src/avx512_intrinsics.rs:160-167** — Unsound AVX-512: store through pointer derived from an immutable local (fast_sin_cos_avx512)
  - Casting `*const T` obtained from `&T`/an immutable place to `*mut T` and writing through it is undefined behavior in Rust: the borrow checker and codegen are allowed to assume the array is never mutated, so the compiler may keep indices_array as all-zeros, making every output sin_out/cos_out read lut[0]. The function is unsound regardless of whether the store currently appears to work. (It is pres
  - *Фикс:* Declare `let mut indices_array = [0i32; 16];` and store via `indices_array.as_mut_ptr() as *mut __m512i`, or use a gather/scalar extraction that does not alias an immutable local.
- **src/sat_if_signal.rs:224-227, 785-788** — L2P code phase reset to 0 every ms (pilot_length forced to 1) → P-code restarts each ms, second half of buffer never played
  - The intent of start_code_phase is to accumulate code phase continuously across ms boundaries (this is exactly the continuity fix described in CLAUDE.md for other signals). With pilot_length=1 the wrap collapses the anchor to 0 each ms, so the L2P P-code repeats with a 1 ms period instead of advancing continuously, and the upper half of the generated 20460-chip P-code buffer is dead code. A receive
  - *Фикс:* Set pilot_len for L2P to the true code length used for wrapping (e.g. data_length = 20460 or the real per-epoch chip count) instead of 1, and wrap start_code_phase by that length so the anchor advances continuously across ms boundaries.


---

## Low (детально)

- **src/ifdatagen.rs:2292** — Redundant always-true modulo masks a latent stale-param bug if PARAM_UPDATE_INTERVAL_MS ever changes
  - Currently harmless because the interval is 1 (params updated every ms, as intended to match C++). But the surrounding code couples this with `full_update = ms_offset == 0` (line 2311): only block-boundary ms get a full update_satellite_params (nav-bit reload + last_nav_bit_index reset). If a future maintainer raises PARAM_UPDATE_INTERVAL_MS to reduce Kepler cost (the comment at 2282 implies that w
  - *Фикс:* Either remove the modulo and call the update unconditionally every ms (documenting that per-ms update is mandatory), or, if decoupling is desired, make push_sat_param_for_ms run every ms and the expensive Kepler get_satellite_param only every interval, and remove the `full_update` block-boundary cou
- **src/satellite_param.rs:409** — get_satellite_cn0 produces NaN CN0 if elevation goes negative in auto-CN0 mode
  - Satellites enter the visible arrays when elevation >= 5° mask at selection time (src/ifdatagen.rs:1811 etc.), but Elevation is recomputed every ms in get_satellite_param during the run. For a long run or a moving receiver a borderline satellite can cross below 0°, yielding NaN here. The resulting CN0=0 then feeds ComputationCache::update -> cached_amp = sqrt(10^(0/10)/fs) (a finite but wrong tiny 
  - *Фикс:* Clamp the argument before sqrt: use `Elevation.sin().max(0.0).sqrt()`, or skip/zero the fade contribution when Elevation <= 0.
- **src/ifdatagen.rs:2229-2231** — Initial AGC gain uses noise variance sigma^2 instead of complex noise power 2*sigma^2
  - The initial AGC gain initial_agc_gain = target_rms / total_rms_with_noise (line 2235-2236) therefore underestimates the true complex RMS by a factor sqrt(2) for the noise-dominated case (with noise_sigma=1.0 and small signal, true block RMS magnitude is ~sqrt(2)=1.414 vs the formula's ~1.0). The first written block is scaled ~41% too hot. The closed-loop RMS controller (line 2423-2436) measures ac
  - *Фикс:* Use the complex noise power: `total_rms_with_noise = (total_rms_amplitude*total_rms_amplitude + 2.0*noise_sigma*noise_sigma).sqrt();` so the initial gain matches the controller's steady-state target.
- **src/ifdatagen.rs:1195, 1206-1208** — select_global_epochs_and_fill ignores the BDS −14 s offset and uses leap=false target
  - Each error (18 s leap on the target, 14 s BDT offset on BDS) shifts that system's absolute time used to score the global candidate. Because candidates are scored by SV count in a ±7200 s window with sum|Δt| tie-break (src/ifdatagen.rs:1210-1284), a ~14-32 s bias is small relative to 7200 s and rarely changes the winner, but the per-system absolute axes are not on a common, exact scale, so the 'uni
  - *Фикс:* Build tgt_abs with leap=true and include the −14 s in bds_abs (bds_abs = (w+1356)*604800 + t - 14) so all three systems and the target share one exact GPS-seconds axis.
- **src/json_interpreter.rs:2691 (consumed by src/satellite_param.rs:660,696)** — GPS mean motion uses EARTH_GM (3.986004418e14) instead of WGS84 GM (3.986005e14)
  - satellite_param.rs uses this `eph.n` directly in the Kepler mean-anomaly propagation (`mk = M0 + (n + alpha/2)*delta_t`, line 660) and in Ek_dot/velocity (line 696). Using μ that is 5.82e7 m³/s² too small biases n by ~Δμ/(2μ)≈7.3e-8 relative, accumulating an along-track orbit/phase error over toe→t. WGS_SQRT_GM exists precisely for the GPS case but is unused.
  - *Фикс:* For the GPS parser use `WGS_SQRT_GM*WGS_SQRT_GM` (3.986005e14) for n; keep EARTH_GM (3.986004418e14) for Galileo, and CGCS2000 μ for BeiDou.
- **src/satellite_param.rs:827 (vs coordinate.rs:426-431); caller get_visible_satellite line 80** — Local gps_sat_pos_speed_eph never reports ephemeris expiry (always returns true)
  - An expired ephemeris (|t-toe|>2h) yields large orbit error yet would be accepted as visible. This silently weakens the visibility gate. Impact is limited because get_visible_satellite has no live caller in the current pipeline (visibility is filtered elsewhere), but the function's contract is violated and any future use is unsafe.
  - *Фикс:* Mirror coordinate.rs: return `delta_t.abs() <= 7200.0` at the end of the local gps_sat_pos_speed_eph (or delete the duplicate and call coordinate::gps_sat_pos_speed_eph).
- **src/satellite_param.rs:822-823 (vs coordinate.rs:420-421)** — BeiDou GEO acceleration Coriolis term uses 2× factor inconsistent with reference
  - The two implementations of the same GEO earth-rotation acceleration compensation disagree by a factor of 2, so at most one matches the C++/ICD reference. The acceleration path is only exercised via get_sat_pos_vel under USE_POSITION_PREDICTION (currently false and with no live callers), so this is presently dead code, but it is a latent correctness divergence if prediction is ever enabled.
  - *Фикс:* Drop the `2.0 *` to match coordinate.rs (single Coriolis factor), or remove the unused prediction/acceleration code path entirely.
- **src/gnsstime.rs:359-400 vs 278-335** — utc_to_glonass_time_corrected (1996 epoch) is not the inverse of glonass_time_to_utc (1992 epoch); incompatible Day/LeapYear encodings
  - The two functions use different Day semantics (absolute-from-epoch vs within-cycle) and different epoch years (1996 vs 1992), so they are not inverses. Verified numerically: utc_to_glonass_time_corrected(2025-06-05) -> Day=10750, LeapYear=7; glonass_time_to_utc of that returns 2023-12-(invalid). The original utc_to_glonass_time (1992) DOES round-trip with glonass_time_to_utc. Impact is currently l
  - *Фикс:* Make the Day encoding consistent: either have utc_to_glonass_time_corrected emit day-within-cycle + N4 (matching glonass_time_to_utc), or provide a corrected inverse glonass_time_to_utc_corrected that decodes the cumulative-day/1996 form. Add a round-trip unit test for the _corrected pair.
- **src/json_interpreter.rs:700-737** — read_nav_file_limited desyncs when a per-system cap is reached: phantom svid=0 ephemerides
  - Verified: read_contents_time on a space-leading data line returns Some(0) (char at index 1 is whitespace -> svid=0, no None early-return). Unlike read_nav_file_filtered, read_nav_file_limited has no dedicated `' '`-only skip arm and actively binds space to the GPS arm. Triggered whenever max_per_system is hit mid-file with mixed-system records remaining (e.g. tests/rinex_parser_baseline.rs uses sm
  - *Фикс:* Add an explicit `' ' => continue` arm before the `'G'` arm (so continuation lines after a skipped record are ignored, matching read_nav_file_filtered), and/or have read_contents_time return None when the SVID position is whitespace instead of Some(0).
- **src/json_interpreter.rs:3198, 3272, 2533, 2553** — Fortran lowercase 'd' exponent not converted, parses to 0.0
  - `f64::from_str("1.234d-09")` errors; the code only replaces "D". The `.unwrap_or(0.0)` then masks the failure, so a clock/orbit term becomes zero with no error. This is a real robustness gap for valid RINEX produced by lowercase-d writers, though the bundled files use lowercase 'e'.
  - *Фикс:* Replace both cases: `line.replace('D', "E").replace('d', "e")` (or uppercase the whole numeric token before parsing).
- **src/json_interpreter.rs:1085, 1546-1553** — skip_ephemeris_lines over-consumes for SBAS (skips 7 lines, SBAS records have 3)
  - Verified against Rinex_Data/BRDC00IGS_R_20251560000_01D_MN.rnx: consecutive S21 records are 4 lines apart (1 header + 3 data). skip_ephemeris_lines loops `0..7`. Within an all-SBAS block this only over-skips harmlessly, but the record immediately after the last SBAS record (if within 7 lines) has its header consumed and is lost.
  - *Фикс:* Skip exactly the number of data lines for the record type (3 for SBAS/GLONASS, 7 for GPS/BDS/GAL), or rely on the `' '`-leading skip behaviour and remove the explicit over-skip for SBAS.
- **src/bcnavbit.rs:185** — MEO reference-axis condition is correct as written; the CLAUDE.md 'documented bug' fix would break it (test/doc gap)
  - The genuine defect is a documentation/test gap, not the line-185 condition. The B-CNAV ephemeris round-trip test (bcnavbit.rs:710) exercises only one MEO satellite (svid 23, flag 3) with tgd_ext all zero, so it neither validates GEO/IGSO satellites, the IODC accumulation path, nor the tgd_b2a source mix-up, and it cannot detect the very issue the documentation claims exists. The misleading note ri
  - *Фикс:* Correct the CLAUDE.md note to reflect the +1 flag mapping (MEO arrives as flag==3, so line 185 is correct). Add round-trip tests covering GEO (flag 1) and IGSO (flag 2) reference axes, non-zero tgd_ext to cover the B2a TGD path, and repeated set_ephemeris calls to cover the data[3] accumulation.
- **src/inavbit.rs:768-769** — Dead code: ephdata[16] &= 0x03ffffff immediately overwritten by assignment
  - The masking on line 768 has no effect because line 769 unconditionally replaces the whole value. This signals confused intent around preserving vs. overwriting word-5 bits (related to the iono-clobber finding) and is at minimum dead code.
  - *Фикс:* Remove line 768, or if preservation of bits [25:0] was intended, change line 769 to `ephdata[16] = (ephdata[16] & 0x03ffffff) | 0x14000000;`.
- **src/memory_code.rs:13** — Galileo E6 memory code array is 16 words short → PRN 50 pilot last 511 chips silently zeroed
  - Verified numerically: E6 file has 15984 hex literals (confirmed via grep/python). 100*160=16000. For svid 50 pilot, slice length 144 < 160; get_memory_sequence chips where word_index (i*32 + (j>>5)) >= 144 — i.e. sector i=4 chips j>=512 — fall into the out-of-bounds `else` and become 0. That zeros 511 of the 5115 pilot chips, corrupting the only channel E6 actually transmits (satellite_signal.rs:5
  - *Фикс:* Pad the E6 memory code file to a full 16000 words (add the missing 16 words of the 100th code) and update the declared size to [u32; 16000], or sourced from the true Galileo E6 memory-code table. Until fixed, restrict Galileo svid to the supported range and bounds-check the slice length before calli
- **tests/prn_generation_baseline.rs:226-261** — Test coverage gap: memory-code, GLONASS, and E5a/E5b codes are only length/balance-checked, missing known-answer and periodicity tests
  - balance_ratio<0.10 is satisfied by a code that merely repeats a balanced 341-chip sub-sequence (G3) and by a code with 511 zeroed chips (E6 svid 50). A correct PRN baseline must assert auto-correlation peak vs sidelobe bounds and full-period non-repetition; the current suite asserts neither, so it cannot detect these correctness failures.
  - *Фикс:* Add tests that (a) generate E6/E1 for the maximum supported svid (50) and assert no all-zero tail and correct length; (b) assert GLONASS G1/G3 and E5a/E5b codes are non-periodic over their full length (no sub-period divides the length) and meet auto-correlation sidelobe bounds; (c) pin at least one 
- **src/navbit.rs:408** — navbit::crc24q_encode panics on empty bit_stream (no guard, unlike crc24q.rs)
  - The two crc24q_encode implementations diverge in robustness for the same documented contract; navbit.rs will panic on `bit_stream[0]` for an empty input while the canonical version returns 0. This is a latent panic, not a wrong CRC value.
  - *Фикс:* Mirror the crc24q.rs guard in navbit.rs, or delegate navbit::crc24q_encode to crate::crc24q::crc24q_encode (as inavbit.rs/fnavbit.rs already do) to keep a single, guarded implementation.
- **tests/pvt_decode_baseline.rs:76-77** — verify_inav_crc24q has no test against real I-NAV encoder output (CRC verifier effectively unverified)
  - The function is documented as 'Verifies CRC24Q for Galileo I-NAV data' and is a public API consumed by external callers, yet a wrong padding/length/word-order assumption (e.g. passing data_len != 196, or a different word layout than encode_data) would still pass the empty-array tests. crc24q_compute uses length_bits.div_ceil(32)*4 so the CRC covers padding bits past the nominal data length; any ca
  - *Фикс:* Add a test that builds an I-NAV page via INavBit::set_ephemeris/get_frame_data (or directly reconstructs encode_data[0..7] as in inavbit.rs:508-515), computes the reference CRC with crate::crc24q::crc24q_encode(&encode_data, 196), and asserts verify_inav_crc24q(&encode_data, 196, reference_crc) == t
- **src/nav_decode.rs:111-113** — extract_bits computes (1u32 << length) which is undefined-shift / returns 0 for length >= 32
  - It is a latent correctness/UB trap: the function is a general-purpose public bit extractor (re-exported and used by tests in tests/pvt_decode_baseline.rs). The current decoders never call it with length >= 32 (they assemble 32-bit fields by OR-ing two sub-extractions and use shifts, max extract length is 24), so no live decode path is wrong today. But the masking idiom is incorrect at the boundary
  - *Фикс:* Make the mask width-safe: let mask = if length >= 32 { u32::MAX } else { (1u32 << length) - 1 }; (word >> position) & mask, mirroring the existing guard in sign_extend. Optionally debug_assert!((1..=32).contains(&length)).
- **src/fastmath.rs:152-184** — fast_gaussian_noise spare path returns imag=0 and only one Gaussian value, corrupting noise statistics
  - A complex noise sample is meant to carry two independent Gaussian samples (I and Q). On the spare path the imaginary (Q) component is forced to 0.0, so every other returned sample has zero quadrature noise. This both wastes the cached sample and produces a non-circular, statistically wrong I/Q noise distribution. The function has no production callers today (only `generate_noise_block` is used, wh
  - *Фикс:* Either remove the spare-caching scheme and always return both fresh values (like `generate_noise_block` does), or have the spare path return a fresh second Gaussian for the imag component instead of 0.0.
- **src/fastmath.rs:51-52, 152-178** — Unsynchronized static mut HAS_SPARE/SPARE in fast_gaussian_noise is a data race if called concurrently
  - The crate uses rayon for satellite-parallel generation. If `fast_gaussian_noise` were called from multiple threads it would be an unsynchronized read-modify-write of shared mutable statics — a data race (UB). It is safe today only because the production noise path uses the thread-local `generate_noise_block` and `fast_gaussian_noise` has no callers. The risk is that the public, seemingly-safe API 
  - *Фикс:* Make the noise generator stateless (as `generate_noise_block` already is) or move the spare cache into thread-local storage / a per-thread RNG struct; do not keep generator state in `static mut`.
- **src/ifdatagen.rs:3166-3177** — LUT-based fast_sin/fast_cos used for satellite orbital position introduces ~km-scale error near elevation mask
  - The module doc (fastmath.rs:15) states LUT trig error <1e-4 'sufficient for IF generation', but here it is reused for geometry that drives the elevation mask. A ~km position error can flip a satellite just above/below the configured elevation mask, changing the visible-satellite set. This is a coarse simplified path (circular M0+n*t, not full Kepler), so the LUT is a secondary error source, but `a
  - *Фикс:* Use standard `mean_anomaly.cos()/.sin()` (and `i0.sin()`) for orbital geometry in `get_visible_satellite`; reserve the LUT for the per-sample carrier loop where the error tradeoff is intended.
- **src/cuda_acceleration.rs:107-127** — GPU C/A PRN kernel is numerically incorrect and O(n^2), only saved by being benchmark-only
  - Real GPS L1 C/A uses a fixed G1 and a G2 tapped at two SV-specific positions; this kernel ignores the SV phase selection entirely, so the generated 'PRN' is identical for every satellite and does not match any real C/A code. It is also O(chip_idx) per sample (quadratic over a code period). The code comment admits it is 'упрощенная для демонстрации', and the kernel is only invoked from benchmarks (
  - *Фикс:* Either document the function clearly as a non-functional benchmark stub, or implement correct per-SV G2 tap selection with incremental LFSR state (not full recomputation per sample).
- **src/avx512_intrinsics.rs:67-94, 137-168, 177-198, 207-234, 242-270** — Dead unused unsafe AVX-512 and lut_sin_cos APIs lack any tests, masking the soundness bugs above
  - Because none of these unsafe functions are exercised, the UB in `fast_sin_cos_avx512` (immutable-array write), the missing `cfg(target_arch)` compile breakage, and the unverified integer-phase LUT indexing in `lut_sin_cos` go undetected. `parallel_satellite_processing_avx512` additionally computes `cos = sqrt(1 - sin^2)` where `sin_approx = new_phases * prn` can exceed 1, making the argument to `_
  - *Фикс:* Add unit tests that validate the SIMD/LUT outputs against scalar references (guarded by `is_x86_feature_detected!`), or delete the unused unsafe functions to remove the latent UB surface entirely.
- **src/almanac.rs:321-358** — read_almanac_galileo dead-guard: missing <issueDate> leaves week=0 instead of triggering the negative-week early return
  - When the header has no parseable <issueDate>, parsing continues with ref_week = 0. get_almanac_galileo then enforces `(ref_week & 3) != data` (line 467) using ref_week=0, which silently rejects (returns None) every almanac whose WNa low-2-bits are non-zero, and accepts wrong ones otherwise; alm.week is also set to 0. The function returns a misleading count and corrupted/empty almanac set instead o
  - *Фикс:* Track whether issueDate was successfully parsed with an explicit bool (or initialize Week to a sentinel like -1) and return 0 early when no reference week was obtained, rather than relying on Week < 0.
- **src/almanac.rs:543-559** — GLONASS almanac: parts[10] parsed into alm.dt then immediately overwritten; dt_dot never populated
  - The parsed column at parts[10] (rate-of-period / dt_dot in the GLONASS almanac record) is computed and thrown away on the very next data-bearing assignment, so it has no effect. Meanwhile dt_dot is left at its Default 0.0, meaning any consumer relying on the draconian-period rate gets no data. This is a silent data-loss / unfinished-parse bug; either the read of parts[10] is wrong or the assignmen
  - *Фикс:* Assign parts[10] to alm.dt_dot (and keep alm.dt = period - 43200.0), or remove the stray line 550 read if parts[10] is not dt_dot. Confirm the column mapping against the almanac format.
- **src/almanac.rs:537** — GLONASS almanac time uses 1992 epoch (utc_to_glonass_time) — the epoch CLAUDE.md flagged as a fixed bug
  - CLAUDE.md documents that the GLONASS epoch must be 1996, not 1992, and that the 1992 reference was a fixed time bug. Here the almanac path still uses the 1992-epoch variant. For the current use (relative newest-wins ordering) the constant epoch offset cancels, so it is not presently fatal, but the stored leap_year/day are on the wrong epoch and any absolute interpretation (or future use of these f
  - *Фикс:* Use utc_to_glonass_time_corrected for almanac timestamping, and make read_almanac_glonass's ordering key consistent with whichever Day convention (day-in-cycle vs cumulative) that function produces.
- **src/almanac.rs:614-626** — get_almanac_from_ephemeris_glonass: unguarded division by eph.vz / pos_vel.vz can produce NaN/Inf
  - If eph.vz or pos_vel.vz is 0.0 (or near zero), the division yields Inf/NaN, which propagates into t and then into glonass_sat_pos_speed_eph on the next iteration, corrupting alm.t, alm.lambda and the whole derived almanac silently (no error returned). Unlike the a<=0 guard at line 654, there is no protection here. Default-constructed or degenerate ephemerides (vz=0) would crash the numerics rather
  - *Фикс:* Guard the divisions (e.g. if pos_vel.vz.abs() < eps return an invalid almanac with flag=0), or clamp/skip the Newton step when |vz| is below a threshold, matching the explicit validity checks already used elsewhere in this function.
- **src/gps_pilot/mod.rs:284-286** — Unbounded carrier-phase accumulation loses precision at nonzero IF (diverges from main generator's bounded model)
  - With the shipped binaries if_freq_hz = 0.0, so carrier_phase only grows with Doppler (~1.5e6 cycles over 300 s) and the f64 ULP there (~2.3e-10 cycles) keeps error sub-millidegree — fine. But GpsPilotConfig/GalileoPilotConfig expose if_freq_hz, and a realistic nonzero IF (e.g. 4.092 MHz) over 300 s drives carrier_phase to ~1.2e9 cycles where the f64 ULP is ~2.4e-7 cycles. Each per-sample increment
  - *Фикс:* Reduce the accumulator to its fractional part each iteration (or each ms) before calling sin_cos, e.g. keep `carrier_phase = carrier_phase.fract()` (or `rem_euclid(1.0)`), or split into an integer-cycle counter plus a fractional remainder as sat_if_signal.rs does. This keeps the sin/cos argument in 
- **src/gps_pilot/mod.rs:232-233** — Receiver ECEF/LLA recomputed every millisecond per satellite inside the hot loop
  - These two trig-heavy coordinate conversions are invariant across the entire generation, yet they are evaluated total_ms × num_satellites times (e.g. 300000 × ~10 = 3M+ redundant lla_to_ecef+ecef_to_lla round-trips for a 300 s GPS run, and proportionally more for the combined gnss_pilot). The receiver is static, so the result is bit-identical every call. This is pure wasted work in the innermost lo
  - *Фикс:* Compute receiver_ecef and receiver_lla once before the chunk/satellite loops (they are already computed once at generate() startup, lines 89-90) and pass the precomputed values into the closures instead of recomputing per ms.
- **src/gnsstime.rs:105 (also satellite_param.rs:249)** — i32 overflow in GPS-week-to-seconds arithmetic (Week * 604800) for weeks > 3550
  - i32::MAX is 2,147,483,647; `Week * 604800` exceeds it once Week > 3550 (3550*604800 = 2.147e9). The multiplication is performed in i32 and only then cast to u32, so for weeks beyond 3550 it overflows: panic in debug builds, silent wraparound in release (overflow-checks not enabled in the release profile in Cargo.toml). A wrapped value yields a wrong leap-second count, corrupting UTC<->GPS conversi
  - *Фикс:* Cast to a wider type before multiplying, e.g. `((gnss_time.Week as i64) * 604800 + (gnss_time.MilliSeconds / 1000) as i64) as u32` (and likewise in satellite_param.rs:249), or compute in u32/i64 throughout.
- **src/fastmath.rs:151-184** — fast_gaussian_noise drops the imaginary component and is not thread-safe (static mut spare)
  - Box-Muller produces two independent Gaussian samples per draw; the function maps them to one complex sample on the fresh branch (real=u1*f, imag=u2*f) but on the spare branch returns only a real value with imag=0, so the cached path produces a complex sample with zero quadrature noise — statistically incorrect complex Gaussian noise. The `static mut` access is also a data race if ever called from 
  - *Фикс:* Generate two real samples per call and return them as (real, imag) without a global spare, or remove the function in favor of generate_noise_block which already uses a local rng. If kept, make it produce a proper complex sample and avoid `static mut`.
- **src/satellite_param.rs:245 (also 1549)** — BeiDou time adjustment subtracts 14000 ms without week-boundary wrap
  - During the first 14 seconds of a (GPS-scale) week, MilliSeconds < 14000, so adjusted_time.MilliSeconds goes negative and satellite_time becomes a small negative seconds-of-week value. delta_t = satellite_time - toe in gps_sat_pos_speed_eph then becomes a large magnitude that may not be corrected by the single ±302400 week-wrap, yielding a wrong satellite position near the week boundary. The genera
  - *Фикс:* After subtracting 14000 ms, if MilliSeconds < 0 add 604800000 and decrement Week (proper week borrow), matching utc_to_bds_time's boundary handling.
- **src/coordinate.rs:1001-1005** — cis_to_cts writes satellite position into the acceleration output instead of acceleration
  - The acceleration array is populated with rotated position values, not acceleration, so any caller requesting GLONASS acceleration via this path gets meaningless data (off by ~7 orders of magnitude and wrong physical quantity). The live GLONASS path passes acc=None (satellite_param.rs:1572), so this is currently dead, but the code is incorrect and a latent trap if acceleration output is ever enable
  - *Фикс:* Either remove the acc output from cis_to_cts (since the propagated state carries no acceleration), or thread the true acceleration from the RK4 integrator (state[6..9] are the lunisolar perturbation accelerations) and rotate those instead of position.
- **src/satellite_param.rs:245** — BeiDou time adjustment can underflow seconds-of-week in get_satellite_param
  - satellite_signal.rs:300-303 explicitly protects the analogous BDS subtraction with `if MilliSeconds < 0 { += 604800000 }`, but the satellite_param.rs path has no such guard. The negative satellite_time produces a delta_t = satellite_time - eph.toe that is off by a full week; coordinate.rs:108 only adds 604800 once when delta_t < -302400, which mostly recovers it, but the edge case is unguarded and
  - *Фикс:* After line 245 add the same wrap guard used in satellite_signal.rs: `if adjusted_time.MilliSeconds < 0 { adjusted_time.MilliSeconds += 604800000; }`.
- **src/ifdatagen.rs:4570-4591** — Dead Galileo ephemeris-copy helper mixes GST-week and GPS-week absolute times
  - Because the 1024-week offset shifts the target far below all candidate epochs, |epoch_abs - target_abs| no longer measures temporal nearness; the min-diff selection degenerates to picking the epoch with the smallest week/toe rather than the one closest to the simulation time. The live epoch selector (select_global_epochs_and_fill, line 1196) correctly treats Galileo week as GPS-aligned, so only th
  - *Фикс:* Use utc_to_gps_time (GPS-aligned) for target_abs in this helper, matching the live path, or delete the dead helper to avoid future reuse of the buggy logic.
- **src/bcnav1bit.rs:658-669** — B-CNAV1 SatType field encoded as 0/1/2 instead of ICD 1/2/3
  - BDS-SIS-ICD-B1C defines the SatType field encoding as 01=GEO, 10=IGSO, 11=MEO (i.e. 1/2/3); 00 is reserved. Emitting 0/1/2 mislabels every satellite type in the generated B-CNAV1 frame and is also inconsistent with the rest of the codebase, which uses flag/sat_type = 1/2/3 (types.rs:413 to_gps_ephemeris maps sat_type+1; bcnavbit.rs:185 checks ==3 for MEO). Only affects receivers decoding the SatTy
  - *Фикс:* Encode sat_type as 1 (GEO), 2 (IGSO), 3 (MEO) to match the ICD and the rest of the codebase.
- **src/ifdatagen.rs:1054, 4099** — Two cur_time initialization paths disagree on leap-second handling
  - The same UTC instant maps to two different GPS times (an 18 s difference in 2025) depending on entry point, and the epoch-selection target is offset by 18 s from the generation clock. For signal correlation an 18 s offset is large (it shifts the satellite geometry/Doppler relative to the selected ephemeris epoch), though the ±7200 s epoch window masks it for epoch selection. The inconsistency is a
  - *Фикс:* Pick one convention (leap seconds true, consistent with load_config and gps_pilot) and use it everywhere cur_time and the epoch-selection target are derived from UTC.
- **src/inavbit.rs:449-451, 497** — Galileo I/NAV WN word-0 write uses unwrapped week and is unmasked
  - If start_time.MilliSeconds ever exceeds one week (>=604800000), start_time.Week is one less than the correct GST week used for tow, so WN and TOW in word 0 disagree. Additionally, with no 12-bit mask, a GST week >= 4096 (or any high bits) would overflow into the TOW bits of data_vec[3]. Word 5 masks correctly, so word 0 is inconsistent. Not currently triggered for normal inputs (MilliSeconds < 1 w
  - *Фикс:* Use start_time_local.Week for the WN field and apply `& 0xfff` as in the word-5 path: data_vec[3] |= ((((start_time_local.Week - 1024) & 0xfff) as u32) << 20) + tow as u32.
- **src/sat_if_signal.rs:398, 793, 531** — Dead AVX-512 / fast-path signal generators contradict documented 'always enabled' acceleration
  - CLAUDE.md states 'AVX-512 fast path (всегда включены)' / 'always enabled', but the AVX path is dead code; all output goes through the scalar f64 sin_cos path. This is a real code/intent mismatch and unmaintained complexity (e.g. the AVX path still uses a separate scalar nav_value and a 1024-entry PrnCache that is incorrect for non-1023-length codes), which is a latent trap if it is ever re-enabled
  - *Фикс:* Either wire the AVX path back into the loop with correct redirects (it already self-redirects non-1023 codes to the cached path) and add a test, or delete the dead functions and the PrnCache to remove the misleading 'enabled' claim.
- **src/gnsstime.rs:236-244** — SubMilliSeconds stores a sub-second (not sub-millisecond) fraction in utc_to_gps_time, double-counting fractional seconds
  - If a start time has a fractional second (e.g. Second=12.3456), 345 ms are added to MilliSeconds (losing 0.6 ms to truncation), yet SubMilliSeconds is set to 0.3456 which downstream is interpreted as 0.3456 ms — the fractional second is both folded into ms AND re-added as a sub-ms term, an inconsistent, double-counted time. The field is documented as sub-millisecond [0,1). Currently dormant because
  - *Фикс:* Set SubMilliSeconds to the sub-millisecond remainder only: `let frac_ms = (utc_time.Second % 1.0) * 1000.0; MilliSeconds += frac_ms.floor(); SubMilliSeconds = frac_ms.fract();` so MilliSeconds + SubMilliSeconds reconstructs the time exactly.
- **src/fastmath.rs:51-52, 152-184** — fast_gaussian_noise uses unsynchronized `static mut HAS_SPARE/SPARE` (data-race-prone, currently dead)
  - If fast_gaussian_noise were ever invoked from the rayon-parallel satellite loop (par_iter_mut in ifdatagen.rs:2358), concurrent read-modify-write of the static mut SPARE/HAS_SPARE would be a data race (UB) and would also corrupt the Box-Muller pairing. It is sound only because it is currently unreachable; the dead unsafe global state is a latent hazard.
  - *Фикс:* Delete fast_gaussian_noise (and the SPARE/HAS_SPARE statics) since generate_noise_block already provides race-free local-rng Box-Muller, or move the spare state into a thread-local / per-call local so it cannot be shared across threads.
- **src/bcnavbit.rs:655-671** — BeiDou MEO reference-axis test only exercises a flag value the data pipeline never produces (test gap)
  - The encoder's correctness for real MEO satellites depends entirely on the to_gps_ephemeris +1 offset matching the bcnavbit `==3` check, but no test asserts that linkage end-to-end; the existing test bypasses it by injecting flag=3 directly, leaving the mapping unguarded.
  - *Фикс:* Add a test that constructs a BeiDouEphemeris with sat_type = BDS_SAT_MEO (=2), runs to_gps_ephemeris(), feeds the result to BCNavBit::set_ephemeris, decodes via decode_bcnav_ephemeris, and asserts the recovered semi-major axis matches the input (catching any drift between the sat_type→flag mapping a
- **tests/orbit_geometry_baseline.rs:199-228 and 131-160** — No test covers relativistic correction or get_glonass_satellite_param scaling
  - Both substantiated correctness bugs (zero relativistic correction; GLONASS double-scaling) pass the existing suite, demonstrating the assertions are too loose / point at the wrong code path. A regression here would ship undetected.
  - *Фикс:* Add a test that compares get_satellite_param's pseudorange/clock to an independent reference (or asserts the relativistic term is nonzero for e>0), and a test that calls get_glonass_satellite_param and asserts Elevation/range are physically plausible (range ~19000-26000 km, not ~2.6e10 m).
- **tests/if_output_smoke.rs:100-109** — Main IF generation pipeline correctness is untested by default (only an 8-sample non-zero smoke test runs)
  - With 8 samples there is not even one full 1 ms code period for GPS L1CA (5000 samples/ms), so the test cannot detect code-phase drift, modulo errors (rem_euclid wrapping), carrier-phase discontinuity, Doppler sign errors, AGC/quantization scale errors, or any nav-bit content defect. The repository's own CLAUDE.md history shows these exact classes of bugs (modulo 0x3FF, code-phase 50 ms jumps, IF f
  - *Фикс:* Add a fast (1-2 s, single-PRN or 2-3 PRN, low sample-rate) default-enabled end-to-end test that generates IQ via IFDataGen and asserts acquisition z-score above threshold and correct code phase for at least one satellite. Keep the 30 s versions ignored, but ensure at least one non-ignored test exerc


---

## Отклонено скептиком (8) — НЕ баги или недоказуемо

- **src/json_interpreter.rs** — RINEX time-window filter comment claims leap seconds are applied but code passes false
  - причина: The claim's central reasoning is a misreading. It asserts GPS/GAL time-window filtering is "inconsistent with BDS filtering" because "the BDS target is derived leap-true via utc_to_bds_time." That is false. is_beidou_ephemeris_within_time_window (src/json_interpreter.rs:1483) builds its target with 
- **src/coordinate.rs** — GLONASS RK4 propagation drops the sub-COARSE_STEP residual, leaving up to ~30 s (~100 km) of orbit motion uncomputed
  - причина: The residual-dropping code in src/coordinate.rs:463-486 is factually real: step_number = (delta_t / COARSE_STEP) as i32 truncates toward zero, the loop runs only the integer number of 30 s RK4 steps, and the residual delta_t - step_number*COARSE_STEP (up to ~30 s) is computed on line 486 but bound t
- **src/gnsstime.rs** — Test gap: time conversions are never exercised with fractional seconds or at leap-second boundaries
  - причина: The claim asserts time-conversion tests are never exercised with fractional seconds or via the _corrected GLONASS pipeline. It looked only at the #[cfg(test)] module in src/gnsstime.rs (lines 522-618) and missed the integration test file tests/time_system_baseline.rs, which directly invalidates two 
- **src/bcnavbit.rs** — clock_param[svid][3] (IODC word) uses |= without prior clear, accumulating stale bits on repeated set_ephemeris calls
  - причина: The literal code observation is true: at src/bcnavbit.rs:320 `data[3] |= COMPOSE_BITS!(eph.iodc as u32, 0, 10);` uses `|=` while the sibling first-writes use `=` (data[0] line 314, data[1] line 316, data[2] line 318). data[3] (= clock_param[svid_idx][3]) has no prior `=` clear inside set_ephemeris, 
- **src/crc24q.rs** — Table-based CRC24Q (src/crc24q.rs) does not match the standard MSB-first CRC-24Q used by Galileo/GPS
  - причина: The claim's central assertions are false and inverted. I compiled and ran the actual code from src/crc24q.rs and compared it against a textbook canonical CRC-24Q (polynomial 0x1864CFB, init 0, MSB-first, no reflection — the RTCM/Galileo standard).

FACT 1 — the table is correct, not "stray". Each CR
- **src/pilotbit.rs** — GPS L5 pilot secondary code stored MSB-first but read LSB-first (bit-reversed output)
  - причина: The claim's conclusion is inverted; the emitted L5 pilot overlay is actually CORRECT. Verified facts:

1. Reader (src/satellite_signal.rs:439-440) reads LSB-first: chip c -> bit (c&0x1f) of word c/32.
2. SECONDARY_CODE_L5 = [0x72b20] (src/pilotbit.rs:46). Read LSB-first the chip sequence is chip0..1
- **src/ifdatagen.rs** — GLONASS satellite index `eph.n as usize - 1` can underflow/panic on RINEX-controlled satellite number
  - причина: The claimed underflow/out-of-bounds panic at src/ifdatagen.rs:3774 and :3823 (`self.glo_sat_param[eph.n as usize - 1]`) is UNREACHABLE because `eph.n` is provably constrained to 1..=24 before it can reach `create_glonass_signals`.

Data-flow proof:
- `create_glonass_signals` only iterates `self.glo_
- **src/json_interpreter.rs** — BeiDou MEO ephemerides in the nav-diagnostics path use the wrong reference semi-major axis (~14760 km encoded orbit error)
  - причина: The claim's mechanism is factually wrong at every step.

(1) Wrong lines: json_interpreter.rs:2701 is inside parse_gps_ephemeris and :3150 is inside parse_galileo_ephemeris. Neither is a BeiDou parser. The BeiDou parser is parse_beidou_ephemeris (line 2945).

(2) The BeiDou parser DOES derive sat_ty