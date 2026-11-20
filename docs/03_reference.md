# 3. Reference

Detailed description of the notebooks, the files they read and write, and the Python modules. For the underlying model see [01_theory.md](01_theory.md); for a walk-through see [02_user_guide.md](02_user_guide.md).

Units: SI unless a name says otherwise (`_um` micrometres, `_in` inches, `_mm_s` mm s⁻¹). Time is seconds since the mixer was switched off. $K$ is in m² s⁻¹.

---

## 3.1 Repository layout

```
README.md
requirements.txt
0_ADV_data_tank_v01.ipynb        template: prepare an ADV decay record
1_OBS_timeseries_v04.ipynb       template, step 1: prepare an OBS record (decay only, t = 0 at mixer-off)
2_SettlingCalibration_v01.ipynb  template, step 2: calibrate the turbulence on a size standard
3_SettlingFit_v01.ipynb          template, step 3: fit one or many runs with a calibration file
tools/settling_tools.py          fitting, calibration, calibration files, ADV shape, PSD handling, settling law
tools/settling_solver.py         numerical model (finite-volume advection-diffusion solver)
tools/UKDE_despike_v02.py        ADV despiking (used by the ADV notebook)
tools/setup_work.py              creates or tops up work/ from work_template/ and the templates
work_template/                   tracked starter files copied into work/: sample OBS and ADV records,
                                 example PSD files, example calibration
work/                            your working folder (ignored by git; created by tools/setup_work.py):
  0_…3_*.ipynb                     your copies of the notebooks
  data_OBS/                        OBS records and step 1 outputs (*_time_series_decayONLY.csv)
  data_ADV/                        ADV records (raw/, filtered/, figures/)
  data_PSD/                        particle size distribution csv files
  calibrations/                    calibration files (written by step 2, read by step 3)
  results/                         step 3 outputs (results csv, fit_plots/)
docs/                            this documentation
```

Everything the notebooks read and write is inside `work/`, which is not part of the repository. The first cell of every notebook (template or copy) finds the project folder and changes into `work/` (creating it with `setup_work` if it does not exist), so **all file paths in the notebook inputs and in this reference are relative to `work/`**: `data_OBS/…` means `work/data_OBS/…`. `python tools/setup_work.py` never overwrites existing files; rerun it after `git pull` to add new sample files or templates.

---

## 3.2 File formats

### OBS decay csv (steps 2 and 3 input; Step 1 output `*_time_series_decayONLY.csv`)

| Column | Meaning | Used by the fit |
|---|---|---|
| `sec` | time since mixer off [s] | yes |
| `NTU` | OBS reading [NTU] | yes |
| `Co_NTU_avg` | initial (mixer-on) concentration [NTU]; first row | yes |
| `Cbackground_NTU_avg` | background reading after settling [NTU]; first row | yes |
| `DateTime`, `DateTime_values`, `second_running`, `min_running`, `C_mg_l`, `Co_mgL_avg`, `Cbackground_mgL_avg` | bookkeeping and mg/L versions | no |

The fit uses $c = (\mathrm{NTU} - C_{bg})/(C_o - C_{bg})$. Other column names can be given with `load_kw` (see `load_decay`).

### ADV decay csv (step 2 input)

Required columns: `decay_time_s` (s since mixer off), `v1_m_s`, `v2_m_s`, `v3_m_s` (velocity components, m s⁻¹). Other columns (`turb_energy`, amplitudes, SNR, correlation) are ignored. Column names can be changed with `t_col` and `vel_cols` in `adv_shape`. The record should reach at least 60 s (for the noise floor) and use a constant sampling rate.

### PSD csv (step 3 input)

One row per size bin: a diameter column and an amount column.

| Option | Default | Meaning |
|---|---|---|
| `d_col`, `p_col` | first and second column | names or 0-based positions |
| `p_type` | `'fraction'` | `'fraction'`: amount per bin (%, volume or mass; rescaled to sum 1). `'cumulative'`: cumulative percent or fraction **finer** than the diameter (differenced) |
| `d_units` | `'um'` | `'um'`, `'mm'` or `'m'` |
| `d_min`, `d_max` | none | keep only bins with $d_{min} \le d \le d_{max}$ (µm) |

Bins with zero amount are removed and the diameters sorted. A PSD finer than `n_max_classes` is merged into log-spaced classes (class diameter = amount-weighted mean). A one-row file is a single size.

### Calibration file (json; written by step 2, read by step 3)

| Key | Content |
|---|---|
| `format`, `version` | `"settling_calibration"`, 2 |
| `created` | timestamp |
| `rpm` | mixer rotation rate |
| `geometry` | `h` (water depth, m), `z_obs` (sensor height above the bed, m) |
| `default_model` | `"adv"` or `"exp"`: the model step 3 uses unless its `turb_model` says otherwise |
| `models` | one entry per calibrated turbulence model, `"adv"` (ADV shape; present only if an ADV record was used) and `"exp"` (exponential; always present), each with the keys below |
| `models.*.K0`, `models.*.Kmin` | calibrated values, m² s⁻¹ |
| `models.*.tau` | s, exponential model only (`null` for `"adv"`) |
| `models.*.free` | parameters that were fitted |
| `models.*.rms`, `models.*.rho1` | rms residual and lag-1 residual autocorrelation of the calibration fit |
| `models.*.sigma_ln` | fractional 1-σ of each fitted parameter (autocorrelation-inflated) |
| `models.*.shape` | `"adv"` only: `t`, `shape` (K/K₀ vs time) and `noise` (noise floor of the detrended fluctuation variance) |
| `models.*.K0_vs_standard_density` | table: `rho_standard`, `K0`, `tau`, `Kmin`, `rms` |
| `standard` | size range, number of classes, density, T, S, `coef`, class weights of the standard |
| `files` | names of the OBS and ADV files used (`adv_file` is `null` without an ADV record) |
| `notes` | free text |

Version 1 files (one model, with `turb_mode`, `K0`, … at the top level) are still read; they hold only the model they were made with.

### Results table (step 3 output, `results/SettlingFit_results.csv`)

| Column | Meaning |
|---|---|
| `run`, `obs_file`, `psd_file` | identification |
| `what`, `units` | fitted quantity (`rhos`/`nf`/`ws50`) and its unit (kg/m3, -, mm/s) |
| `value`, `sigma` | fitted value and its 1-σ fit uncertainty |
| `value_K0_low`, `value_K0_high`, `K0_factor` | fitted value with K₀ divided and multiplied by `K0_factor` (present when `do_sensitivity` is true) |
| `w50_mm_s`, `w50_lo_mm_s`, `w50_hi_mm_s` | settling velocity at the median size and the range implied by the 1-σ (for `ws50`: ± its 1-σ) |
| `w_mean_mm_s` | amount-weighted mean velocity |
| `d50_um` | median diameter of the PSD used |
| `rms`, `rho1` | rms residual and lag-1 residual autocorrelation |
| `turb_model` | turbulence model used, `adv` or `exp` |
| `K0`, `Kmin`, `tau` | calibration values used (`tau` empty for `adv`) |

---

## 3.3 `1_OBS_timeseries_v04.ipynb` (step 1)

**Purpose.** Turn a raw OBS logger file into a time series and a decay-only csv with t = 0 at mixer-off.

**Sections.** Find the OBS files; read in the data; time corrections; plot the reformatted data; extra processing for decay analysis (averaging windows, mixer-off time); trim the time series; optional conversion to mg/L.

**Input variables** (set in several cells, top to bottom).

| Variable | Default | Meaning |
|---|---|---|
| `dir_data`, `ext` | `'data_OBS/'`, `'*.csv'` | folder and pattern of the raw files |
| `j` | `filenames.index('MicroSphere_75_90_SAMPLE')` | the file to process: by name, or an index from the printed list |
| `hr_cor`, `min_cor`, `sec_cor` | 0, −1, 2 | added to the logger's first timestamp (the logger time of a row is that of its last sample) |
| `Co_start_s`, `Co_end_s` | – | averaging window for $C_o$ (s from the start of the record) |
| `Cbackground_start_s`, `Cbackground_end_s` | – | averaging window for $C_{background}$ (s) |
| `mixer_on` | – | when the mixer was switched on (s); only used to define `mixer_off` |
| `mixer_off` | `mixer_on + 4*60` | when the mixer was switched off (s); sets t = 0 |
| `m`, `b` | – | optional NTU → mg/L calibration, $C_{mass} = m\,C_{NTU} + b$ |

**Raw file.** Logger csv (TOA5): four header lines, then one row per minute: timestamp, record number, one NTU value per sample. The notebook moves the third column to the end before stacking the samples.

**Outputs** (in `dir_data`).

* `<name>_time_series.csv`: `DateTime`, `DateTime_values`, `second_running` (s from the start of the record), `min_running`, `NTU`, and `C_mg_l` if the mg/L cell is run.
* `<name>_time_series_decayONLY.csv`: the same columns for t > `mixer_off`, plus `sec`, `Co_NTU_avg`, `Cbackground_NTU_avg`, and the mg/L versions if converted. Format in section 3.2.

---

## 3.4 `2_SettlingCalibration_v01.ipynb` (step 2)

**Purpose.** Calibrate both turbulence models on a size standard, $K_0$ and $K_{min}$ for the ADV shape and $K_0$, $\tau$ and $K_{min}$ for the exponential, and write them to one calibration file. Without an ADV record (`adv_file = None`) only the exponential model is calibrated.

**Sections.** 1 input; 2 standard decay record; 3 ADV shape (skipped without an ADV record); 4 calibrate each model; 5 K₀ vs the density of the standard, per model; 6 save; 7 optional second-standard test, per model.

**Input cell.**

| Variable | Default | Meaning |
|---|---|---|
| `standard_obs_file` | – | OBS decay csv of the standard |
| `adv_file` | – | ADV decay csv, same mixer rate; `None` if there is no ADV record (only the exponential model is calibrated) |
| `calibration_file` | `calibrations/cal_60rpm.json` | output file |
| `rpm` | 60 | mixer rate; stored and checked in step 3 |
| `notes` | `''` | free text stored in the file |
| `water_depth_in`, `sensor_above_bottom_in` | 7.0, 2.5 | tank geometry in inches; the sensor height is measured **up** from the bottom |
| `T_C`, `S_ppt` | 19.4, 0.0 | water temperature (°C) and salinity (ppt) |
| `coef` | 5 | settling-law coefficients (1–5, section 3.8) |
| `d_lo_um`, `d_hi_um`, `n_classes` | 75, 90, 7 | size range of the standard and the number of equal-width size bins |
| `rho_standard` | 1200 | density of the standard (kg m⁻³) |
| `psd_weights` | `None` | `None`: uniform in diameter; or `n_classes` weights from a measured PSD |
| `default_model` | `'adv'` | model step 3 uses by default; becomes `'exp'` when `adv_file = None` |
| `free_adv` | `('K0','Kmin')` | parameters fitted for the ADV model; `K0` and/or `Kmin` |
| `free_exp` | `('K0','tau','Kmin')` | parameters fitted for the exponential model; any subset of `K0`, `tau`, `Kmin` |
| `K0_start`, `tau_start`, `Kmin_start` | 1e-4, 15, 5e-6 | starting values (and fixed values for parameters not fitted) |
| `adv_window_s` | 2.0 | length of each ADV averaging window (s) |
| `adv_noise` | `None` | ADV noise floor (m² s⁻²); `None` estimates it from the late record |
| `adv_noise_from_s` | 60 | the noise floor is estimated from the record after this time (s) |
| `adv_cut_s` | 30 | the shape is set to zero after this time (s); $K_{min}$ carries the weak mixing beyond it |
| `adv_normalize` | `'first'` | `'first'`: shape = 1 in the first window; `'max'`: shape = 1 at its maximum |
| `density_offsets` | (−50, −25, 0, 25, 50) | density offsets for the K₀ table; `()` skips it |
| `validate`, `val_obs_file`, `val_d_lo_um`, `val_d_hi_um`, `val_rho` | `False`, … | optional prediction of a second standard |

**Outputs.** ADV shape plot; fit and residual plot with each model; a K₀-vs-density table per model; the calibration file; optionally the second-standard prediction and its rms for each model. Run time about three minutes with both models (most of it the density tables).

---

## 3.5 `3_SettlingFit_v01.ipynb` (step 3)

**Purpose.** Fit density, fractal dimension or the median-size settling velocity for any number of runs using a calibration file.

**Sections.** 1 input; 2 calibration (loads the model(s) used, checks geometry and rate, plots $K(t)$ of each); 3 fit all runs (fit, plot, save one plot per run); 4 results table (saved as csv); 5 class velocities of one run.

**Input cell.**

| Variable | Default | Meaning |
|---|---|---|
| `calibration_file` | `calibrations/cal_60rpm.json` | from step 2 |
| `rpm` | 60 | rate of these runs, checked against the file (`None` skips the check) |
| `runs` | example | list of dicts: `name`, `obs_file`, `psd_file`, plus optional per-run overrides of any default below |
| `turb_model` | `None` | turbulence model from the calibration file: `'adv'`, `'exp'`, or `None` for the file's `default_model` |
| `what_to_fit` | `'rhos'` | `'rhos'`, `'nf'` or `'ws50'` |
| `T_C`, `S_ppt` | 19.4, 0.0 | water temperature and salinity |
| `coef` | 5 | settling-law coefficients (1–5, section 3.8) |
| `nf` | 3.0 | fractal dimension when fitting density (3 = solid) |
| `dp_um` | 8.0 | primary particle size (µm); no effect for $n_f = 3$ |
| `rhos_for_nf` | 2650 | density used when fitting $n_f$ |
| `beta` | 1.3 | $d$-exponent of the velocity spectrum for `'ws50'` |
| `gain_exp` | 0.0 | OBS signal per unit amount ∝ $d^{-\mathrm{gain\_exp}}$ |
| `psd` | dict | PSD csv layout (`d_col`, `p_col`, `p_type`, `d_units`, `d_min`, `d_max`) |
| `n_max_classes` | 12 | merge a finer PSD into at most this many classes |
| `t_fit_start`, `t_fit_end` | 0.0, `None` | fit window (s since mixer off) |
| `do_sensitivity` | `True` | also fit with K₀ divided and multiplied by `K0_factor` (about 3× longer) |
| `K0_factor` | 1.4 | K₀ uncertainty factor |
| `geometry_override` | `None` | a `settling_tools.Geometry` to replace the geometry stored in the calibration |
| `results_csv`, `plot_dir` | `SettlingFit_results.csv`, `fit_plots` | output file and plot folder names, inside `results_folder` |
| `results_folder` | `'results/'` | folder for all step 3 outputs; created if missing |

Per-run overrides recognised: `turb_model`, `what_to_fit`, `T_C`, `S_ppt`, `coef`, `nf`, `dp_um`, `rhos_for_nf`, `beta`, `gain_exp`, `psd`, `t_fit_start`, `t_fit_end`, `K0_factor`.

**Outputs.** The results table (section 3.2), one plot per run in `results_folder + plot_dir` and in the notebook, and a class-velocity table for one run. Run time about 25 s per run with sensitivity, about 10 s without.

---

## 3.6 `tools/settling_tools.py`

### Constants and tables

* `DIFF_STEP = 1e-3`: relative step (in the log-parameters) for the finite-difference Jacobians of all least-squares fits.
* `COEFS`: settling-law coefficient sets $(b_1, b_2)$: 1 → (20, 1.258) flocs matching Stokes at the primary particles; 2 → (120, 0) viscous flocs; 3 → (18, 0.866) sieve-sized sand; 4 → (20, 0.91) nominal-diameter natural sand; 5 → (18, 0.548) smooth sphere.
* Parameter bounds in the fits: settling velocity 1e-5 to 5e-2 m s⁻¹; $K_0$ 1e-7 to 1e-1 m² s⁻¹; $\tau$ 0.5 to 900 s (1 to 900 s in `fit_ws` with `fit_turbulence`); $K_{min}$ 5e-9 to 1e-2 m² s⁻¹; excess density ($\rho_s - \rho_w$) 1 to 6000 kg m⁻³; $n_f$ 1.2 to 3.0. Bounds for the calibration can be changed with `bounds={'K0': (lo, hi), ...}`.

### Water and settling law

| Function | Description |
|---|---|
| `density(To, Sppt)` | water density (kg m⁻³) from temperature (°C) and salinity (ppt) |
| `visc(To)` | kinematic viscosity (m² s⁻¹) |
| `settlingvel(d_um, dp_um, nf, coef, rhos, T, S)` | settling velocity (m s⁻¹); `d_um` may be an array |
| `invert_ws(ws, d_um, which, T, S, coef, rhos=None, nf=None, dp_um=8.0)` | `which='rhos'` returns the density (needs `nf`); `which='nf'` returns the fractal dimension (needs `rhos`, `dp_um`); NaN if there is no solution |
| `class_ws(d_um, which_known, T, S, coef, rhos=None, nf=3.0, dp_um=8.0)` | velocities of an array of diameters for a known density (`which_known` is unused) |

### Data and geometry

| Function / class | Description |
|---|---|
| `load_decay(path, t_col="sec", ntu_col="NTU", co_col="Co_NTU_avg", cb_col="Cbackground_NTU_avg")` | returns dict `t`, `ntu`, `Co`, `Cb`, `c` |
| `t50_estimate(t, c, window=7)` | time the smoothed normalised concentration first falls below 0.5; a starting guess only (residual turbulence makes it overestimate the velocity) |
| `Geometry(h, z_obs)` | dataclass: `h` water depth, `z_obs` sensor height above the bed (m); properties `z_s = h − z_obs` and `fall = h − z_obs` |

### Size distributions

| Function | Description |
|---|---|
| `size_classes(d_lo, d_hi, n=7, weights=None)` | `n` equal-width bins over a size range; returns diameters (µm) and fractions |
| `load_psd(path, d_col=None, p_col=None, p_type="fraction", d_units="um", d_min=None, d_max=None)` | read a PSD csv; returns diameters (µm) and fractions summing to 1 |
| `rebin_psd(d_um, p, n_max=12)` | merge into at most `n_max` log-spaced classes |
| `signal_weights(d_um, p, gain_exp=0.0)` | OBS signal fraction per class, $p_i (d_i/d_{50})^{-\mathrm{gain\_exp}}$ renormalised |

### Turbulence from the ADV

| Function | Description |
|---|---|
| `adv_shape(path, t_col="decay_time_s", vel_cols=("v1_m_s","v2_m_s","v3_m_s"), window_s=2.0, noise=None, noise_from_s=60.0, t_cut_s=30.0, smooth=3, normalize="first")` | returns dict `t`, `shape`, `noise`, `fluct` (variance per window), `umax`. `noise` is in m² s⁻² of the **detrended** fluctuation variance (not the raw sum of squares of the velocities, which includes mean offsets). `normalize="max"` scales to the maximum instead of the first window |
| `shape_K(shape, K0, Kmin=0.0)` | callable $K(t) = K_0\,\mathrm{shape}(t) + K_{min}$, shape held at its first value before the first window and zero after the last |

### Forward model

| Function | Description |
|---|---|
| `model_c(t, ws, K0, tau, Kmin, geom, N=600, Kfun=None)` | normalised sensor concentration, one class. `Kfun` (callable) replaces the exponential $K(t)$ |
| `model_c_classes(t, ws, p, K0, tau, Kmin, geom, N=300, Kfun=None)` | several classes (velocities `ws`, signal fractions `p`) in one integration |
| `model_c_psd(t, fit, geom, K0, tau, Kmin, shape=None, N=250)` | model curve for the classes and velocities of a `fit_psd` or `fit_ws_classes` result |
| `predict_standard(t, geom, cal, d_lo, d_hi, rho, T, S, coef, n_classes=7, weights=None, N=250)` | model decay of a standard with a stored calibration, no fitting |

### Fitting and calibration

| Function | Description |
|---|---|
| `fit_ws(t, c, geom, K0, tau, Kmin=5e-8, ws0=None, t_start=0.0, t_end=None, fit_turbulence=False, N=400, shape=None)` | fit the velocity of a single class. With `shape` the turbulence is $K_0\,\mathrm{shape} + K_{min}$ and `fit_turbulence` fits $K_0$ only. Returns dict `ws`, `ws_sigma`, `K0`, `tau`, `rms`, `rho1`, `n`, `t`, `c`, `resid`, `success` |
| `fit_psd(t, c, geom, d_um, p, which, T, S, coef, K0, tau=30.0, Kmin=5e-8, shape=None, nf=3.0, dp_um=8.0, rhos=2650.0, gain_exp=0.0, t_start=0.0, t_end=None, N=250)` | fit density (`which='rhos'`) or fractal dimension (`'nf'`) for a size distribution. Returns dict `which`, `value`, `sigma`, `ws` (class velocities), `d_um`, `p`, `weights`, `d50`, `w50`, `w50_range`, `w_mean`, `rms`, `rho1`, `resid`, `t`, `c`, `success` |
| `fit_ws_classes(t, c, geom, d_um, p, K0, tau, Kmin=5e-8, ws50_0=None, t_start=0.0, t_end=None, beta=2.0, N=250, shape=None)` | fit the median-size velocity with $w_{s,i} = w_{s,50}(d_i/d_{50})^\beta$. Returns dict `ws50`, `ws50_sigma`, `ws`, `d50`, `rms`, `rho1`, `resid`, `t`, `c`, `success` |
| `fit_run(obs_file, psd_file, geom, cal, T, S, coef, what="rhos", nf=3.0, dp_um=8.0, rhos_for_nf=2650.0, beta=1.3, gain_exp=0.0, psd_kw=None, n_max_classes=12, t_start=0.0, t_end=None, K0_factor=1.4, do_sensitivity=True, N=250, name=None, load_kw=None)` | one run from files to results: returns `(row, fit, data)`: the results-table row, the fit dict and `dict(t, c, d_um, p, name)`. `psd_kw` are the `load_psd` options, `load_kw` the `load_decay` options |
| `sensitivity_table(t, c, geom, K0_list, tau_list, d_um, T, S, coef, nf=3.0, dp_um=8.0, Kmin=5e-8, ws0=None, t_start=0.0, t_end=None, N=200, shape=None)` | refit the single-class velocity for each $K_0$ (and $\tau$); a DataFrame |
| `sensitivity_psd(t, c, geom, d_um, p, which, T, S, coef, K0_list, tau=30.0, Kmin=5e-8, shape=None, **kw)` | as `fit_psd` for a list of $K_0$ values; a DataFrame |
| `calibrate_turbulence(t, c, geom, ws_known, free=("K0","tau"), K0=1e-4, tau=10.0, Kmin=5e-8, t_start=0.0, t_end=None, N=300, bounds=None)` | calibrate the exponential turbulence on a single-size standard of known velocity |
| `calibrate_turbulence_classes(t, c, geom, ws, p, free=("K0","tau"), K0=1e-4, tau=10.0, Kmin=5e-8, t_start=0.0, t_end=None, N=250, bounds=None, shape=None)` | calibrate on a size-range standard (class velocities `ws` and fractions `p` fixed). With `shape`, `free` may hold only `'K0'` and/or `'Kmin'`. Returns dict `K0`, `tau`, `Kmin`, `rms`, `rho1`, `free`, `resid`, `t`, `c`, `sigma_ln` (fractional 1-σ of each free parameter) |
| `calibrate_models(t, c, geom, ws, p, shape=None, free=None, K0=1e-4, tau=15.0, Kmin=5e-6, N=250)` | calibrate every model the inputs allow: `'exp'` always, `'adv'` when `shape` is given. `free` is `{model: parameters}`, defaulting to `DEFAULT_FREE` (`'adv'`: K0, Kmin; `'exp'`: K0, tau, Kmin). Returns `{model: calibrate_turbulence_classes result}` with `turb_mode` added |
| `cal_Kfun(cal, shape=None)` | callable $K(t)$ of one calibrated model (a `calibrate_models` entry or a `load_calibration` result) |
| `model_c_cal(t, ws, p, cal, geom, shape=None, N=250)` | model decay of classes `ws`, `p` with one calibrated model |
| `calibration_sensitivity(t, c, geom, d_lo, d_hi, n_classes, rhos_list, T, S, coef, shape=None, free=("K0","Kmin"), K0=1e-4, tau=15.0, Kmin=5e-6, weights=None, N=200)` | recalibrate for each assumed density of the standard; DataFrame `rho_standard`, `K0`, `tau`, `Kmin`, `rms` |

### Calibration files

| Function | Description |
|---|---|
| `save_calibration(path, cals, geom, shape=None, rpm=None, standard=None, sensitivity=None, files=None, notes="", default_model=None)` | write the json file (creates the folder). `cals` is `{model: result}` from `calibrate_models` (a single `calibrate_turbulence_classes` result is also accepted), `shape` the `adv_shape` dict (required for `'adv'`), `sensitivity` `{model: DataFrame}`; `default_model` defaults to `'adv'` when calibrated, else `'exp'` |
| `load_calibration(path, model=None)` | read one model: `'adv'`, `'exp'`, or `None` for the file's default. Dict with `K0`, `Kmin`, `tau`, `tau_used` (30 s placeholder for `'adv'`, where it is unused), `shape` (arrays, or `None`), `turb_mode` (the model), `available_models`, `geom` (a `Geometry`), `rpm` and the stored diagnostics. `ValueError` if the model is not in the file. Reads version 1 files |
| `check_calibration(cal, geom=None, rpm=None, tol=0.02)` | list of warning strings: depth or sensor height differing from the file by more than `tol` (2 %), or a different mixer rate; empty if consistent |

### Fit uncertainty

All least-squares fits report the Jacobian-based 1-σ scaled by the residual variance and inflated by $(1+\rho_1)/(1-\rho_1)$ with $\rho_1$ the lag-1 autocorrelation of the residuals (clipped to 0–0.95). It is a fit uncertainty only.

---

## 3.7 `tools/settling_solver.py`

| Function | Description |
|---|---|
| `solve_single(w, Kfun, h, z_s, t_eval, N=400, D=1.0, K0_init=None, sigma_s=0.0, c_init=None, rtol=1e-8, atol=1e-10, return_profiles=False, scheme="limited")` | one class with depth-mean initial concentration 1. `Kfun(t)` is the diffusivity; `z_s` is the sensor **depth below the surface**; `D` = 1 for a perfectly depositional bed (0 reflecting); `K0_init` (default `Kfun(0)`) sets the stirred initial profile; `sigma_s` is a Gaussian sensor footprint (0 = one cell); `scheme="limited"` (second-order Van Leer) or `"sg"` (first-order exponential fitting); `return_profiles` also returns the grid, profiles and deposited fraction |
| `solve_classes(ws, p, Kfun, h, z_s, t_eval, N=300, D=1.0, K0_init=None, sigma_s=0.0, rtol=1e-7, atol=1e-9, k=None, return_classes=False)` | several classes in one integration; returns $\sum_i k_i c_i(z_s,t)$ (`k` = per-class sensor gain, default 1) and per-class signals if requested |
| `solve_mixture(ws, p, k, Kfun, h, z_s, t_eval, **kw)` | same as `solve_classes` by looping over `solve_single` (reference) |
| `K_exp(K0, tau, Kmin=0.0)` | callable exponential $K(t)$ |
| `theta_exp(K0, tau, Kmin=0.0)` | callable $\int_0^t K\,dt'$ for the exponential model |
| `stirred_profile(w, K0, h, N, mean=1.0)` | stirred steady initial profile |
| `sensor_weights(zc, z_s, sigma_s, dz)` | sensor weighting over the cells |
| `erfc_front(w, theta, z_s, t)` | far-boundary analytic approximation $\tfrac12\,\mathrm{erfc}[(wt - z_s)/\sqrt{4\theta}]$ for checks |

`z` is measured downward from the free surface; the Python-level fitting functions take a `Geometry` (sensor height above the bed) and convert.

---

## 3.8 Settling-law coefficient (`coef`) guidance

| `coef` | Intended for |
|---|---|
| 1 | flocs: matches Stokes at the primary particles; use with `nf` < 3 |
| 2 | viscous flocs |
| 3 | sieve-sized sand |
| 4 | natural sand by nominal diameter |
| 5 | smooth spheres (also the size standards) |

---

## 3.9 Reference results (75–90 µm standard, 1200 kg m⁻³, 60 rpm, fresh water at 19.4 °C)

Benchmarks from the sample data (`MicroSphere_75_90_SAMPLE.csv`, `ADV_02_60rpm_SAMPLE.dat`). Use them to check an installation and as a guide to what normal values look like for this tank.

**Turbulence models compared** (calibration on the standard; the last two rows are the two models step 2 stores by default):

| Turbulence model | K₀ (m² s⁻¹) | τ (s) | K_min (m² s⁻¹) | rms |
|---|---|---|---|---|
| exponential, K₀ and τ free | 9.5e-5 | 34 | 5e-8 (fixed) | 0.0814 |
| exponential, K₀, τ and K_min free (default) | 1.05e-4 | 18 | 4.0e-6 | 0.0799 |
| ADV shape, K₀ free | 1.3e-4 | — | 5e-8 (fixed) | 0.090 |
| ADV shape, K₀ and K_min free (default) | 1.2e-4 | — | 6.2e-6 | 0.081 |

The models fit about equally well (rms 0.080–0.090 against sensor noise of about 0.04), so the OBS decay alone does not discriminate between them; this is why the ADV shape is used. Freeing K_min improves the rms from about 0.090 to 0.081, although K_min itself is only weakly determined (fractional 1-σ about 0.2).

**K₀ against the assumed density of the standard** (ADV shape, K₀ and K_min free): 1150 kg m⁻³ → 5.0e-5, 1175 → 7.9e-5, 1200 → 1.2e-4, 1225 → 1.9e-4, 1250 → 3.9e-4 m² s⁻¹, with rms 0.080–0.084 throughout. K₀ changes by a factor of about 2.5–3 for ±50 kg m⁻³.

**Sensitivity of the fit to K₀:** halving or doubling K₀ moves the fitted velocity of the standard from 0.60 to 0.76 mm s⁻¹ (implied density 1171–1217 kg m⁻³ against a true 1200).

**Step 3 on the standard** (also the calibration run, so this only confirms the code path): density 1198 ± 3 kg m⁻³, w₅₀ = 0.698 mm s⁻¹, rms 0.081; with K₀ divided and multiplied by 1.4 the density is 1186 and 1209 kg m⁻³. A one-row PSD of 82.5 µm gives 1197 kg m⁻³. With the exponential model (`turb_model = 'exp'`) the 75–90 µm fit gives 1200 ± 3.5 kg m⁻³, w₅₀ = 0.703 mm s⁻¹, rms 0.080. A 7-class fit takes about 5–10 s, about 25 s with the two K₀ sensitivity fits.

**Other mixer rates:** the 30 rpm ADV record gives the same K₀ and K_min when its shape is used, although it has about 7× less initial energy and a broader plateau before decaying. Use the ADV record that matches the mixer rate of the OBS runs.

---

## 3.10 Version notes

* Calibration file version 2 holds both turbulence models (ADV shape and exponential), or only the exponential when there is no ADV record; step 3 chooses with `turb_model`. Version 1 files still load.
* v01 of the two notebooks introduces the calibration file, the multi-run fit and the repository layout. `tools/settling_tools.py` was previously `settling_single.py`.
* Earlier single-size notebooks (`SettlingFit_SingleClass_v01/v02`, `SettlingFit_v03`, `SettlingFit_PSD_v01`) are superseded by steps 2 and 3 and are not part of the repository.
