# Audit: `PGA_PGV_PGD_Fourier_v3.m` (PGV / PGD and baseline correction)

Files in this folder:

| File | Purpose |
|---|---|
| `PGA_PGV_PGD_Fourier_v3.m` | original script (unchanged, for reference) |
| `process_ground_motion.m` | corrected processing function |
| `PGA_PGV_PGD_Fourier_v4.m` | driver script: read, process, print PGA/PGV/PGD, fc sweep, plots, export |
| `test_process_ground_motion.m` | tests against records whose exact answer is known |

## How it was checked

I built a synthetic record whose exact acceleration, velocity and displacement are known
(a 0.3 Hz Gaussian-windowed pulse: PGV = 56.55 cm/s, PGD = 27.90 cm). I then added the
usual baseline errors: a 0.3 cm/s² offset, a 0.8 cm/s² baseline step at t = 18 s and
noise. I ran the v3 methods and the corrected pipeline on it with Python/SciPy and with
Octave 8.4 plus the signal package. MATLAB itself was not available.

| Method (v3) | PGV [cm/s] | PGD [cm] | Comment |
|---|---|---|---|
| **exact** | **56.55** | **27.90** | |
| M1 raw cumtrapz | 52.6 | 549.8 | quadratic drift |
| M2 zero mean | 62.8 | 138.1 | the step still causes drift |
| M3 detrend | 58.0 | 46.1 | +65 % |
| M4 filter + detrend v, d (fc = 0.1) | 56.3 | 27.9 | OK, but see findings 3 and 4 |
| M4 when filter is off (`fc = 0` branch) | 54.3 | 32.0 | +15 %, and v(0) = 2.1 cm/s |
| **v4 `process_ground_motion`** | **56.5** | **28.06** | v(0) ≈ 0, a/v/d consistent |

With the 0.1 Hz filter on, v3 Method 4 gives about the same PGD as the corrected pipeline.
**The large PGD differences you see come mainly from findings 1, 2 and 5 below (units,
filter settings and the cutoff frequency), not from the integration itself.**

## Findings (most important first)

### 1. Units and labels are wrong. This is the most likely cause of PGD mismatches
* The comment says the acceleration is in `cm2/s`, the plot labels say `m/s²`, and the
  velocity label says `cm/s`.
* `displacement1_mm = displacement1;` converts nothing. If the acceleration is in cm/s²,
  the displacement is in **cm**, so every value labelled mm is **10× off**.
* PEER NGA files (`.AT2`) store acceleration in **g**. If your file is in g, every
  velocity and displacement is **980.665× too small** until it is converted.
* v4 has an explicit `opt.units` setting ('g', 'cm/s2', 'm/s2'). It always outputs
  cm/s², cm/s and cm, and has a unit test.

### 2. PGA, PGV and PGD are never computed
The script is named PGA_PGV_PGD but never computes `max(abs(...))`. v4 returns
`R.PGA`, `R.PGV` and `R.PGD` together with the times at which they occur.

### 3. `detrend(velocity)` and `detrend(displacement)` in Method 4 are not physical
* `detrend` also removes the mean. Afterwards v(0) ≠ 0 and d(0) ≠ 0, although the ground
  starts at rest.
* After those two lines, the velocity is no longer the integral of the acceleration you
  plot, and the displacement is no longer the integral of that velocity. The three
  series you report disagree with each other.
* It hides drift instead of fixing it. It also removes any genuine permanent
  displacement (fling step). Several near-fault Kocaeli stations recorded fling.
* When the filter is on, the damage is small (under 1 % on the test record). When the
  filter is off, PGD is 15 % too high.
* **Fix:** apply every correction to the **acceleration** only, then integrate once.
  v4 does this. If a residual displacement drift remains, v4 has an optional
  PEER-style correction (`opt.disp_poly_order`). It fits t²…tᴺ to the displacement
  and subtracts the polynomial's **second derivative from the acceleration**, so a, v
  and d stay consistent.

### 4. The Butterworth filter uses the unstable `[b,a]` form
At a normalised corner of 0.1/100 = 0.001, the transfer-function form behaves like this
(checked against the SOS form):

| order | error of `[b,a]` vs SOS (fc = 0.1 Hz, Fs = 200 Hz) |
|---|---|
| 4 | 1.6e-7, acceptable |
| 6 | 1.8e+6, **unstable** |
| 8 | 1e+133, **unstable** |

If your software lets users choose the filter order, `[b,a]` will produce garbage from
order 6 upwards. v4 uses `butter → zpk → zp2sos → sosfilt`.

### 5. PGD depends mostly on the cutoff frequency and the filter type
Test record containing 0.15, 0.3 and 1 Hz components (true PGV = 51.4 cm/s,
PGD = 28.4 cm):

| fc_low [Hz] | acausal PGV / PGD | causal PGV / PGD |
|---|---|---|
| 0.03 | 50.7 / 25.8 | 50.8 / 34.9 |
| 0.05 | 51.1 / 26.9 | 46.4 / 31.0 |
| 0.10 | 49.8 / 25.1 | 39.7 / 25.4 |
| 0.20 | 42.1 / 14.1 | 32.7 / 18.7 |

* A fixed `fc = 0.1` Hz is not safe for a Mw 7.4 near-fault record. Real long-period
  motion lies below 0.1 Hz, so choose fc for each record (from the signal-to-noise
  ratio or the FAS) and **report how PGD changes with fc**. v4 prints this sweep.
* `filtfilt` with `butter(4)` is **effectively order 8, −6 dB at fc**, not order 4,
  −3 dB. When you compare with SeismoSignal, use the **same** filter type, order, corner
  frequencies, and causal or acausal setting. Otherwise PGV and PGD will differ by
  10–30 %, as the table shows. v4 exposes `opt.order`, `opt.causal`, `opt.fc_low` and
  `opt.fc_high` for this reason.

### 6. The filter is run without zero padding
MATLAB's `filtfilt` pads only 3·(nfilt−1) = 12 samples (0.06 s). The transient of a
0.1 Hz filter lasts tens of seconds. Boore (2005) and Boore & Bommer (2005) recommend:
* zero pads with a total length of 1.5·order/fc, split between the two ends;
* integrating the padded record and removing the pads **afterwards**.

v4 follows both. On the test records the effect was small. It matters when strong
motion is close to the start or end of the record.

### 7. `T = 0.005` is hard-coded, but `dt` is read from the file
`Fs` (used by the filter) comes from `T`, while the integration uses `time`. If the file
has another dt, the filter corners are wrong. v4 takes dt from the file and stops with
an error if the sampling is not uniform.

### 8. Fourier spectrum
* `P2(1:L/2+1)` fails for an **odd** number of samples. v4 uses `floor(L/2)+1` and has
  a test for it.
* `|FFT|/L` is a single-sided sinusoid-amplitude spectrum. The Fourier amplitude
  spectrum used in seismology (and, as far as I know, in SeismoSignal) is
  **|FFT|·dt**, in cm/s. v4 returns `R.FAS` this way.

### 9. Minor issues
* Integrating displacement with `cumtrapz` is not exact for piecewise-linear
  acceleration. v4 uses the exact formula d₊ = d + v·dt + (2a + a₊)dt²/6 (option
  `'trapz'` keeps the old behaviour). The difference is below 0.1 %.
* Detrending the whole record before filtering lets the strong-motion part bias the
  fitted trend. v4 can remove the pre-event mean (`opt.pre_event_s`) and fit a
  polynomial of order 0–3, like SeismoSignal's baseline correction.
* `results_header` is never written, and `save -ascii` writes no header. v4 writes a
  labelled file, `results_corrected.txt`.

## Running

```matlab
% MATLAB (Signal Processing Toolbox)   |  Octave:  pkg load signal
test_process_ground_motion   % must print ALL TESTS PASSED
PGA_PGV_PGD_Fourier_v4        % edit the settings block first (file, units, fc)
```

Octave 8.4 output on the synthetic record: PGV 56.48 cm/s and PGD 28.06 cm
(exact values: 56.55 and 27.90). All 8 tests pass.
