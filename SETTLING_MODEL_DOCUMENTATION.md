# Settling-velocity model for a stirred mixing tank — documentation

## Contents of the toolbox

**The three-step workflow**

| Step | Notebook | Purpose |
|---|---|---|
| 1 | `1_OBS_timeseries_v04.ipynb` (yours) | Import a new OBS record, find C_o and C_background, trim to the decay after the mixer stops, reset t = 0. Output: `*_time_series_decayONLY.csv` (columns used downstream: `sec`, `NTU`, `Co_NTU_avg`, `Cbackground_NTU_avg`). |
| 2 | `2_SettlingCalibration_v01.ipynb` | Once per mixer rotation rate: calibrate K₀ and K_min on a size standard of known size range and density, with the ADV decay shape. Output: a calibration file (`work/calibrations/cal_60rpm.json`). |
| 3 | `3_SettlingFit_v01.ipynb` | For any number of runs: fit density, fractal dimension or the velocity at the median size, from an OBS decay csv and a PSD csv, using the calibration file. Output: results csv and one plot per run. |

Most users do step 1, then step 3.

**Code:** `tools/settling_solver.py` (numerical model), `tools/settling_tools.py` (fitting, calibration and calibration files, ADV shape, PSD handling, settling law; renamed from `settling_single.py`). Python 3 with numpy, scipy, pandas, matplotlib.

**Superseded notebooks** (replaced by steps 2 and 3): `SettlingFit_v03.ipynb` and `SettlingFit_PSD_v01.ipynb` (updated to import `settling_tools`, so they still run) and `SettlingFit_SingleClass_v01/v02.ipynb` (import the old `settling_single.py`, which is no longer maintained).

---

## 1. What the code does

You measure the total suspended-sediment signal of an optical backscatter sensor (OBS) in a mixing tank. The mixer runs until the suspension is steady, then is switched off (t = 0). Particles settle past the sensor and the signal decays. Residual turbulence decays over the first 20–35 s and spreads the particles diffusively while they settle.

The code:

1. simulates the decay of concentration at the sensor with a 1-D vertical advection–diffusion model that includes a time-decaying turbulent diffusivity K(t) and the tank boundaries;
2. fits the **settling velocity** w to the measured decay;
3. converts w into **particle density** (solid particles, fractal dimension fixed at 3) or **fractal dimension** (flocs, density fixed) through the Ferguson & Church settling law adapted for fractal aggregates;
4. calibrates the turbulence (K₀, and τ or K_min) on a **size standard** of known size range and density, using the ADV to supply the shape of the turbulence decay.

Inputs you always have: particle size (or size range), the OBS concentration time series, water temperature and salinity, tank geometry. Optionally an ADV record of the same mixer setting.

---

## 2. Model

### 2.1 Governing equation

z is measured **downward from the free surface** (z = 0 surface, z = h bed):

    ∂c/∂t + ∂/∂z ( w c − K(t) ∂c/∂z ) = 0

* w — settling velocity of the particle class (m/s, positive downward).
* K(t) — turbulent diffusivity, uniform in depth, decaying with time (section 2.3).
* Surface (z = 0): zero particle flux.
* Bed (z = h): depositional flux F = D·w·c with D = 1 (particles that reach the bed stay there; D = 0 would be a reflecting bed).

### 2.2 Initial condition

With the mixer running the suspension is in steady state under the stirred diffusivity K₀ = K(0): zero net flux gives

    c(z) ∝ exp( w (z − h) / K₀ )

normalised to the depth-mean. This is almost uniform when w h / K₀ ≪ 1 and bottom-heavy otherwise. For several classes each class has its own profile with its own mean (its mass fraction p_i).

### 2.3 Turbulence K(t)

Two options, selected by what you pass to the functions:

| Option | Form | Parameters |
|---|---|---|
| exponential | K(t) = K_min + K₀ exp(−t/τ) | K₀, τ, K_min |
| ADV shape | K(t) = K_min + K₀ · shape(t) | K₀, K_min; shape(t) from the ADV, shape(0) = 1 |

K₀ is the stirred-state diffusivity. K_min is a small persistent floor (the tank is never truly quiescent); it is an effective model parameter and is far below what the ADV can resolve.

**ADV shape** (`adv_shape`). The ADV point velocities cannot be converted into a tank-scale eddy diffusivity without a length scale, so only the time *shape* of the decay is taken from the ADV; the amplitude K₀ comes from the calibration. Steps: windows of 2 s with 50 % overlap; in each window the linear trend of every velocity component is removed (this strips the mean flow and offsets that the raw sum of squares contains); the fluctuation variance v₁′²+v₂′²+v₃′² is computed; the noise floor (median of the record after 60 s) is subtracted; K is taken proportional to the rms fluctuation velocity, shape ∝ √(variance − noise); a 3-window running mean is applied; the shape is set to 0 after `t_cut_s` (default 30 s, where the signal reaches background); and it is normalised to 1 in the first window.

The assumed scaling K ∝ √E (energy decays with time constant τ_E, so K decays with 2τ_E) is an assumption, not a measurement.

### 2.4 Observation and normalisation

The model returns the concentration in the cell at the sensor depth z_s = h − z_obs, where z_obs is the sensor height above the bottom (a Gaussian footprint is available in the solver, off by default). The sensor signal is normalised like the data:

    c = (NTU − C_background) / (C_o − C_background)

and the model is divided by its own value at t = 0. C_o is the initial plateau and C_background is what the OBS reads when settling has finished (first-row values of `Co_NTU_avg` and `Cbackground_NTU_avg`).

### 2.5 Several size classes

The problem is linear and the classes do not interact, so the sensor signal is the sum over classes of the class concentration at the sensor:

    c_obs(t) = Σ_i c_i(z_s, t),   with each class starting from its stirred profile with mean p_i.

All classes are integrated together in one call (`solve_classes`). A size range such as 75–90 µm is split into equal-width bins (`size_classes`), by default uniform in diameter; pass `weights` for a measured PSD. The OBS is assumed to respond equally to all classes (gain k_i = 1 unless `k` is given to the solver).

### 2.6 Numerics

Finite-volume grid of N cells in depth (N = 200–600 depending on function). Advection uses second-order Van Leer limited face values (little numerical diffusion once K has decayed and the settling front is nearly sharp); diffusion is central. Time integration is scipy BDF (rtol 1e-5 … 1e-7). For several classes a coloured finite-difference Jacobian (5 colour groups) is used, which is about 3× faster than scipy's generic sparse Jacobian. A mass-balance variable tracks deposited mass. The multi-class result matches the sum of single-class runs to about 6e-8.

---

## 3. Settling law and water properties

`settlingvel(d_um, dp_um, nf, coef, rhos, T, S)` returns w in m/s:

    w = g R d^(nf−1) / ( C₁ ν d_p^(nf−3) + C₂ √( g R d^nf d_p^(nf−3) ) ),   R = (ρ_s − ρ_w)/ρ_w

* d — particle (floc) diameter, d_p — primary particle diameter, n_f — fractal dimension (3 = solid sphere-like particle).
* ρ_w from `density(T, S)` (temperature °C, salinity ppt); ν from `visc(T)`.
* `coef` selects (C₁, C₂): 1 flocs (20, 1.258); 2 viscous flocs (120, 0); 3 sieve-sized sand / smooth spheres (18, 0.866); 4 natural sand (20, 0.91); 5 sphere (18, 0.548).

`invert_ws(ws, d_um, which, T, S, coef, rhos=None, nf=None, dp_um=8.0)` solves the law for the unknown: `which='rhos'` gives the density (needs `nf`; solid particles nf = 3), `which='nf'` gives the fractal dimension (needs `rhos` and `dp_um`). It returns NaN if there is no solution in range.

---

## 4. Workflow

### Step 1 — OBS record (`1_OBS_timeseries_v04.ipynb`)
Produces the decay csv. `load_decay(path)` reads `sec` (time since mixer off), `NTU`, `Co_NTU_avg` and `Cbackground_NTU_avg` (first-row values) and returns `t`, `ntu`, `Co`, `Cb` and the normalised concentration `c = (NTU − Cb)/(Co − Cb)`. The other columns of the csv are ignored.

### Step 2 — calibrate the turbulence (`2_SettlingCalibration_v01.ipynb`)
Inputs: the standard's OBS decay csv, the ADV decay csv of the same mixer rotation rate, the standard's size range and density, tank geometry, water, rotation rate. The class velocities of the standard are fixed by its size range and density; K₀ and K_min (ADV mode) or K₀ and τ (exponential mode) are fitted.

The notebook (a) builds the ADV shape, (b) fits and plots the standard, (c) recalibrates for assumed densities of the standard ±50 kg/m³ (the realistic uncertainty of K₀), (d) writes the calibration file, and (e) optionally predicts a **second standard** (different size or density) with the stored calibration and no refitting, which is the real test of the calibration.

### Calibration file (json)
`save_calibration` / `load_calibration`. Contents: format and version, creation time, mixer rpm, tank geometry (water depth, sensor height above the bottom), turbulence model (`adv` or `exp`), K₀, K_min, τ (exponential mode only), which parameters were fitted, rms and lag-1 residual autocorrelation, 1-σ of each fitted parameter, the ADV shape (time and normalised K/K₀ arrays and the noise floor), the standard (size range, density, T, S, settling-law coefficient), the table of K₀ against the assumed density of the standard, the source file names and free-text notes. Step 3 needs only this file, not the ADV csv.

`check_calibration(cal, geom, rpm)` returns warnings when the water depth or sensor height differs from the calibration by more than 2 % or when the rotation rate differs; step 3 prints them.

### Step 3 — fit runs (`3_SettlingFit_v01.ipynb`)
Inputs: the calibration file, then a list of runs. Each run gives a name, an OBS decay csv and a PSD csv; any default (water temperature and salinity, settling-law coefficient, what to fit, fit window, PSD layout) can be overridden per run. Geometry and K(t) come from the calibration file. For each run the notebook fits, repeats the fit with K₀ divided and multiplied by `K0_factor` (default 1.4), plots data, model and residual, and writes one row of the results table (`results/SettlingFit_results.csv`): fitted value and 1-σ, velocity at d₅₀ and its range, mean velocity, d₅₀, rms, lag-1 autocorrelation, K₀, K_min and the fitted values at the two K₀ extremes. See section 4b for the PSD file and fit options.

### Use the same setup throughout
The calibration applies to one mixer rotation rate in one tank setup (depth, sensor position). Water temperature and salinity may differ between runs.

### Minimal script

```python
import sys
sys.path.insert(0, '../tools')   # run from work/, so data paths match the notebooks; the modules live in tools/
import numpy as np, settling_tools as ss

geom = ss.Geometry(h=7*0.0254, z_obs=2.5*0.0254)         # depth and sensor height above bottom [m]
T, S, coef = 19.4, 8.0, 3

std = ss.load_decay('data_OBS/standard_time_series_decayONLY.csv');  run = ss.load_decay('data_OBS/unknown_time_series_decayONLY.csv')
shape = ss.adv_shape('data_ADV/filtered/ADV_02_60rpm_SAMPLE_filtered_decay_only.csv')                # turbulence shape (60 rpm record)

# 1. calibrate on the 75-90 um standard of known density
d, p = ss.size_classes(75, 90, 7)
ws_std = ss.class_ws(d, None, T, S, coef, rhos=1200.0)
cal = ss.calibrate_turbulence_classes(std['t'], std['c'], geom, ws_std, p,
                                      free=('K0', 'Kmin'), K0=1e-4, Kmin=5e-6, shape=shape)
ss.save_calibration('calibrations/cal_60rpm.json', cal, geom, shape=shape, rpm=60)   # step 2 output
cal = ss.load_calibration('calibrations/cal_60rpm.json')                              # step 3 starts here

# 2. fit the unknown with that turbulence
t50 = ss.t50_estimate(run['t'], run['c'])
fit = ss.fit_ws(run['t'], run['c'], geom, cal['K0'], 30.0, cal['Kmin'],
                ws0=geom.fall / t50, shape=cal['shape'])

# 3. back out density (solid) or fractal dimension (flocs)
rho = ss.invert_ws(fit['ws'], d_um=120.0, which='rhos', T=T, S=S, coef=coef, nf=3.0)
nf  = ss.invert_ws(fit['ws'], d_um=120.0, which='nf',   T=T, S=S, coef=1, rhos=2650.0, dp_um=8.0)
print(fit['ws']*1e3, '+/-', fit['ws_sigma']*1e3, 'mm/s')
```

---

## 4b. PSD file and fit options (step 3)

**PSD csv.** One row per size bin, a diameter column and an amount column (mass or volume, any units; rescaled to sum 1). The defaults read the first two columns; set `d_col` / `p_col` for others. `p_type='cumulative'` reads "percent finer" tables, `d_units` accepts um, mm or m, and `d_min` / `d_max` trim the range. A PSD finer than `n_max_classes` (default 12) is merged into log-spaced classes to keep fits fast. A **single size is a PSD file with one row**. Examples are in `data_PSD/PSD_example_75_90.csv`, `data_PSD/PSD_example_75_90_cumulative.csv` and `data_PSD/PSD_example_single_82_5um.csv` (synthetic, to show the layout).

**What is fitted** (`what_to_fit`):

| value | unknown | class velocities |
|---|---|---|
| `'rhos'` | density of solid particles (n_f fixed, default 3) | settling law at each class diameter |
| `'nf'` | fractal dimension of flocs (density `rhos_for_nf` and `dp_um` fixed) | settling law at each class diameter |
| `'ws50'` | velocity at the median size, no density or n_f assumption | w_i = w₅₀ (d_i/d₅₀)^β, β set by `beta` |

The PSD amounts are taken as the fraction of the initial OBS signal in each class; if the sensor responds differently to different sizes set `gain_exp` (signal per unit amount ∝ d^−gain_exp).

```python
cal = ss.load_calibration('calibrations/cal_60rpm.json')
for w in ss.check_calibration(cal, cal['geom'], rpm=60): print(w)
row, fit, data = ss.fit_run('data_OBS/run07_time_series_decayONLY.csv', 'data_PSD/run07_PSD.csv', cal['geom'], cal,
                            T=19.4, S=8.0, coef=3, what='rhos', K0_factor=1.4)
print(row['value'], row['sigma'], row['w50_mm_s'], row['value_K0_low'], row['value_K0_high'])
```

Check on the 75–90 µm standard (it is also the calibration run, so this only confirms the code path): density 1199 ± 3 kg/m³ (known 1200), w₅₀ = 0.662 mm/s, rms 0.080; with K₀ divided and multiplied by 1.4 the density is 1185 and 1211 kg/m³. A one-row PSD of 82.5 µm gives 1198 kg/m³. A run with a 7-class fit takes about 5–10 s; with the two K₀ sensitivity fits about 25 s.

---

## 5. API reference

### 5.1 `tools/settling_tools.py`

**Data and geometry**

* `load_decay(path, t_col="sec", ntu_col="NTU", co_col="Co_NTU_avg", cb_col="Cbackground_NTU_avg")` → dict(`t`, `ntu`, `Co`, `Cb`, `c`).
* `t50_estimate(t, c, window=7)` — first time the smoothed normalised concentration falls below 0.5. Starting guess only; residual turbulence makes it overestimate w.
* `Geometry(h, z_obs)` — dataclass; properties `z_s = h − z_obs` (sensor depth below surface) and `fall = h − z_obs` (distance a particle must fall to clear the sensor).

**Turbulence from the ADV**

* `adv_shape(path, t_col="decay_time_s", vel_cols=("v1_m_s","v2_m_s","v3_m_s"), window_s=2.0, noise=None, noise_from_s=60.0, t_cut_s=30.0, smooth=3, normalize="first")` → dict(`t`, `shape`, `noise`, `fluct`, `umax`). `noise` is in m²/s² of the *detrended* fluctuation variance (not the raw "background energy" of the sum of squares, which includes mean offsets). `normalize="max"` scales to the maximum instead of the first window (use if the record starts with a weak spin-up).
* `shape_K(shape, K0, Kmin=0.0)` → callable K(t).

**Forward model**

* `model_c(t, ws, K0, tau, Kmin, geom, N=600, Kfun=None)` — normalised sensor concentration for one class. `Kfun` overrides K₀, τ, K_min.
* `model_c_classes(t, ws, p, K0, tau, Kmin, geom, N=300, Kfun=None)` — several classes (velocities `ws`, fractions `p`).

**Fitting**

* `fit_ws(t, c, geom, K0, tau, Kmin=5e-8, ws0=None, t_start=0, t_end=None, fit_turbulence=False, N=400, shape=None)` — fits w (and optionally K₀, τ; with `shape` only K₀). Returns `ws`, `ws_sigma`, `K0`, `tau`, `rms`, `rho1`, `resid`, `success`.
* `fit_ws_classes(t, c, geom, d_um, p, K0, tau, Kmin, ..., beta=2.0, shape=None)` — fits one velocity scale for a size range, assuming w_i = w₅₀ (d_i/d₅₀)^β. Not yet in the notebook.
* `calibrate_turbulence(t, c, geom, ws_known, free=("K0","tau"), ...)` — single-size standard, exponential K(t).
* `calibrate_turbulence_classes(t, c, geom, ws, p, free=("K0","tau"), K0, tau, Kmin, t_start, t_end, N=250, bounds=None, shape=None)` — size-range standard; with `shape`, `free` may contain only `'K0'` and/or `'Kmin'`. Returns the fitted values, `rms`, `rho1`, `resid` and `sigma_ln` (fractional 1-σ of each free parameter).
* `sensitivity_table(t, c, geom, K0_list, tau_list, d_um, T, S, coef, ..., shape=None)` — refits w for each K₀ (and τ) and returns a DataFrame with w, implied density and rms. With `shape`, τ is ignored.

**Size distribution (PSD)**

* `load_psd(path, d_col=None, p_col=None, p_type="fraction", d_units="um", d_min=None, d_max=None)` → (`d_um`, `p`), sorted, `p` sums to 1.
* `rebin_psd(d_um, p, n_max=12)` — merge into at most `n_max` log-spaced classes.
* `signal_weights(d_um, p, gain_exp=0.0)` — OBS signal fraction per class.
* `fit_psd(t, c, geom, d_um, p, which, T, S, coef, K0, tau=30, Kmin=5e-8, shape=None, nf=3, dp_um=8, rhos=2650, gain_exp=0, t_start=0, t_end=None, N=250)` — `which='rhos'` or `'nf'`. Returns `value`, `sigma`, `ws` (class velocities), `d50`, `w50`, `w50_range`, `w_mean`, `rms`, `rho1`, `resid`, `weights`.
* `fit_ws_classes(...)` (above) — velocity-scale fit; returns `ws50`, `ws50_sigma`, `ws`, `d50`, `rms`, `rho1`, `resid`.
* `model_c_psd(t, fit, geom, K0, tau, Kmin, shape=None)` — model curve for a fit result.
* `sensitivity_psd(t, c, geom, d_um, p, which, T, S, coef, K0_list, tau=30, Kmin=5e-8, shape=None, **kw)` — refit for each K₀.

**Calibration file and runs**

* `calibration_sensitivity(t, c, geom, d_lo, d_hi, n_classes, rhos_list, T, S, coef, shape=None, free=("K0","Kmin"), K0, tau, Kmin, weights=None, N=200)` — recalibrates for each assumed density of the standard; returns a DataFrame (`rho_standard`, `K0`, `tau`, `Kmin`, `rms`).
* `save_calibration(path, cal, geom, shape=None, rpm=None, standard=None, sensitivity=None, files=None, notes="")` — writes the json file (creates the folder).
* `load_calibration(path)` → dict with `K0`, `Kmin`, `tau`, `tau_used` (30 s placeholder in ADV mode, where it is unused), `shape` (or None), `geom` (a `Geometry`), `rpm`, `turb_mode` and the stored diagnostics.
* `check_calibration(cal, geom=None, rpm=None, tol=0.02)` → list of warning strings (empty = consistent).
* `predict_standard(t, geom, cal, d_lo, d_hi, rho, T, S, coef, n_classes=7, weights=None, N=250)` — model decay of a standard with a stored calibration, no fitting.
* `fit_run(obs_file, psd_file, geom, cal, T, S, coef, what="rhos", nf=3.0, dp_um=8.0, rhos_for_nf=2650.0, beta=1.3, gain_exp=0.0, psd_kw=None, n_max_classes=12, t_start=0.0, t_end=None, K0_factor=1.4, do_sensitivity=True, N=250, name=None, load_kw=None)` → `(row, fit, data)`: one results-table row, the fit dict and the data used. This is what step 3 calls for each run.

**Size classes and settling law**

* `size_classes(d_lo, d_hi, n=7, weights=None)` → bin-centre diameters (µm), mass fractions.
* `class_ws(d_um, which_known, T, S, coef, rhos=None, nf=3.0, dp_um=8.0)` → class velocities (m/s) for a known density.
* `settlingvel`, `invert_ws`, `density`, `visc` — see section 3.

**Constants**: `DIFF_STEP = 1e-3`, the relative step for the least-squares finite-difference Jacobians (the ODE output carries about 1e-5 integration noise, so the scipy default step gives unreliable uncertainties).

### 5.2 `tools/settling_solver.py`

* `solve_single(w, Kfun, h, z_s, t_eval, N=400, D=1.0, K0_init=None, sigma_s=0.0, ...)` — one class, depth-mean initial concentration 1. `K0_init` defaults to `Kfun(0)`.
* `solve_classes(ws, p, Kfun, h, z_s, t_eval, N=300, D=1.0, K0_init=None, sigma_s=0.0, k=None, return_classes=False)` — several classes in one integration; returns the summed signal (and per-class signals if requested).
* `solve_mixture(ws, p, k, Kfun, h, z_s, t_eval)` — same by looping over single-class runs (reference).
* `K_exp(K0, tau, Kmin)`, `theta_exp(...)` — exponential turbulence helpers; `erfc_front(w, theta, z_s, t)` — far-boundary analytic approximation (error O(l_d/z_s)), useful for checking.
* `stirred_profile`, `sensor_weights` — initial profile and sensor footprint.

Note that `z_s` in the solver is the depth below the free surface, but the Python-level fitting functions take `Geometry(h, z_obs)` and convert it.

---

## 6. Uncertainty and what limits the answer

* **± on w** from `fit_ws` is the least-squares uncertainty inflated for lag-1 autocorrelation of the residuals. It does not include K₀, K_min, geometry, size or temperature.
* **K₀ dominates.** In the 75–90 µm standard, halving or doubling K₀ moves the fitted w from 0.56 to 0.74 mm/s (implied density 1169–1221 kg/m³ against a true 1200). K₀ in turn trades off against the assumed density of the standard, about a factor of 2 for ±50 kg/m³. Run the sensitivity table every time.
* **K₀ and τ trade off** when both are fitted from the OBS alone; the ADV shape removes that ambiguity in `'adv'` mode.
* **K_min** is only weakly determined (fractional 1-σ about 0.2) and is far below the ADV noise (7–12 mm/s per sample). It improves the rms from about 0.087 to 0.080.
* **Residual structure**: rms ≈ 0.08 against sensor noise ≈ 0.04, lag-1 autocorrelation ≈ 0.4. Something systematic remains; the standard's actual density or size distribution is the main suspect.
* **Circularity**: if the "unknown" run is the standard itself, the fit just recovers the standard's density (1198 vs 1200 kg/m³ in the example).

## 7. Reference results (75–90 µm standard, 1200 kg/m³, 60 rpm)

| Turbulence model | K₀ (m²/s) | τ (s) | K_min (m²/s) | rms |
|---|---|---|---|---|
| exponential, K₀ and τ free | 8.0e-5 | 32 | 5e-8 (fixed) | 0.0815 |
| exponential, K₀, τ and K_min free | 9.2e-5 | 14 | 3.5e-6 | 0.0798 |
| ADV shape, K₀ free | 1.0e-4 | — | 5e-8 (fixed) | 0.087 |
| ADV shape, K₀ and K_min free | 1.0e-4 | — | 4.4e-6 | 0.080 |

The models fit about equally well (rms 0.080–0.087 against sensor noise of about 0.04); the OBS decay alone does not discriminate between them, which is why the ADV shape is used. K₀ against the assumed density of the standard (ADV shape, K₀ and K_min free): 1150 kg/m³ → 4.2e-5, 1175 → 6.6e-5, 1200 → 1.0e-4, 1225 → 1.5e-4, 1250 → 2.6e-4 m²/s, rms 0.080–0.084 throughout.

The 30 rpm ADV record gives the same K₀ and K_min when its shape is used. It has about 7× less initial energy and a broader plateau before decaying. Use the ADV file that matches the mixer speed of the OBS run.

## 8. Limitations

* 1-D and depth-uniform K; the spin-down circulation of the tank is not represented.
* OBS assumed to respond linearly and equally to all size classes; the sensor footprint is a single cell by default.
* K ∝ √(excess energy) is an assumption; the ADV measures at a point, not the tank average.
* The PSD enters as signal fractions of the initial reading; size-dependent sensor response must be supplied through `gain_exp`. Flocs: the PSD must be the floc size distribution at the time of the run.
* Timing: a fit takes seconds (single class, N = 400) to about ten seconds (7–12 classes, N = 200–250); the calibration notebook takes about 2 minutes with the density table, the fit notebook about 25 s per run with the K₀ sensitivity.

## 9. Troubleshooting

* `WARNING: … differs from the calibration` (step 3) — the water depth, sensor height or mixer rate does not match the calibration; recalibrate for the new setup or correct the input.
* `ValueError: … is not a settling calibration file` — `calibration_file` points to the wrong file.
* `ValueError: ADV record too short…` — the record is shorter than `noise_from_s` (60 s); lower it or pass `noise`.
* `ValueError: with an ADV shape the free parameters can only be 'K0' and/or 'Kmin'` — remove `'tau'` from `free`.
* Fit hits a bound or returns an odd w — give a better `ws0` (use `geom.fall / t50_estimate(...)`, then expect the answer 30–40 % lower).
* Density is NaN — the fitted w is outside the range the settling law can produce for that size (check `d_um`, `coef`, `nf`).
* Results shift when you change `N` — compare N = 200 and 400; differences should be well below the fit uncertainty.
