# Tank settling-velocity toolbox

Estimate the settling velocity of suspended particles from the **total-concentration decay** recorded by an optical backscatter sensor (OBS) in a stirred mixing tank, and convert it to particle density (solid particles) or fractal dimension (flocs).

After the mixer is switched off, particles settle past the sensor while residual turbulence decays over the first 10–35 s. A one-dimensional advection–diffusion model with a decaying turbulent diffusivity $K(t)$ and the tank's free surface and bed describes both effects. $K(t)$ is calibrated once per mixer rotation rate on a **size standard** of known size and density, using an acoustic Doppler velocimeter (ADV) record for the shape of the turbulence decay; settling velocities of unknown populations are then fitted from their decay alone, for any initial **particle size distribution**.

## Workflow

| Step | Notebook | Purpose | How often |
|---|---|---|---|
| 1 | `1_OBS_timeseries_v04.ipynb` | Prepare an OBS record: find $C_o$ and background, keep only the decay after the mixer stops, reset t = 0. Output: `*_time_series_decayONLY.csv`. | every new record |
| – | `0_ADV_data_tank_v01.ipynb` | Prepare an ADV decay record. | with each calibration |
| 2 | `2_SettlingCalibration_v01.ipynb` | Calibrate the turbulence models (ADV shape and exponential, or the exponential alone without an ADV record) on a size standard. Output: a calibration file in `work/calibrations/`. | **once per mixer rate** and tank setup |
| 3 | `3_SettlingFit_v01.ipynb` | Fit density, fractal dimension or the median-size settling velocity for any number of runs. Output: results csv and one plot per run. | every set of runs |

Most users run step 1 and then step 3 with an existing calibration file.

## Quick start

```bash
git clone https://github.com/FluidSedDynamics/TankSettling.git
cd TankSettling
pip install -r requirements.txt
python tools/setup_work.py     # creates work/: copies of the notebooks plus the sample data
jupyter lab
```

Work in the notebooks **inside `work/`**. That folder holds your copies of the notebooks together with all data, calibrations and results, and every path in the notebooks is relative to it. It is ignored by git, so `git pull` updates the code cleanly and nothing you process is ever committed. (Opening a template notebook in the project folder also works in `work/`, and creates `work/` on first run if needed.)

To try it on the included sample (a 75–90 µm standard of known density): run `work/1_OBS_timeseries_v04.ipynb` to make its decay csv, then `work/3_SettlingFit_v01.ipynb`. No inputs need changing; the fit takes about a minute and returns a density of about 1198 kg m⁻³ against a known 1200. See the [user guide](docs/02_user_guide.md#22-installation) for the full sample walk-through.

For your own data, put the records in `work/data_OBS/`, `work/data_ADV/raw/` and `work/data_PSD/`, and change the inputs in the `work/` notebooks. Keep as many copies of a notebook as you like (for example one step 1 notebook per record). `git pull` updates `tools/` for your notebooks too; updated template notebooks reach `work/` only as new files when you rerun `python tools/setup_work.py` (it never overwrites). See `work/README.md` (copied from `work_template/README.md`).

If you change code and push it, also run `pip install nbstripout` and `nbstripout --install` once in your clone, so notebook outputs are never committed.

## Repository layout

```
0_ADV_data_tank_v01.ipynb        ADV preparation (template)
1_OBS_timeseries_v04.ipynb       step 1 (template)
2_SettlingCalibration_v01.ipynb  step 2 (template)
3_SettlingFit_v01.ipynb          step 3 (template)
tools/settling_tools.py          fitting, calibration, calibration files, ADV shape, PSD handling, settling law
tools/settling_solver.py         finite-volume advection-diffusion solver
tools/UKDE_despike_v02.py        ADV despiking (used by the ADV notebook)
tools/setup_work.py              creates or tops up work/
work_template/                   tracked starter files copied into work/: sample OBS and ADV records,
                                 example PSD files, example calibration
work/                            your working folder (ignored by git; created by tools/setup_work.py):
  0_…3_*.ipynb                     your copies of the notebooks
  data_OBS/                        OBS records and step 1 outputs (*_time_series_decayONLY.csv)
  data_ADV/                        ADV records (raw/, filtered/, figures/)
  data_PSD/                        particle size distribution csv files
  calibrations/                    calibration files (written by step 2, read by step 3)
  results/                         step 3 outputs (results csv, fit_plots/)
docs/                            documentation
```

## Documentation

* [Theory of the fit](docs/01_theory.md): the model, the turbulence calibration, estimation and uncertainty, assumptions.
* [User guide](docs/02_user_guide.md): installation, data preparation, a step-by-step walk-through of steps 1 to 3, quality checks, troubleshooting.
* [Reference](docs/03_reference.md): file formats, every notebook input, every function and option, and reference results for the sample standard.

## Things to know before using results

* A calibration is valid for **one mixer rotation rate in one tank setup** (water depth, sensor position). The notebooks warn if these differ.
* The fit's ± is a fit uncertainty only. The turbulence calibration is the largest uncertainty: $K_0$ changes by a factor of about 2.5–3 for ±50 kg m⁻³ in the assumed density of the standard. Step 3 reports the fitted value with $K_0$ divided and multiplied by an uncertainty factor next to every result.
* The model is one-dimensional and does not represent the tank's spin-down circulation. Residuals for the supplied standard (rms about 0.08 against sensor noise of about 0.04) show that something systematic remains.
* Testing the calibration on a second standard of different size or density (an optional section of step 2) is the recommended check.

## Status

Tested with Python 3.13 (numpy 2.5, scipy 1.18, pandas 3.0, matplotlib 3.11).

License and citation: to be added by the repository owner.

## License

MIT; see [LICENSE](LICENSE).
