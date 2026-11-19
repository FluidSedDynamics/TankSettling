# work/ — your working folder

This folder is where you process data. It was created from `work_template/` by `python tools/setup_work.py` (or automatically by the first notebook you ran). Everything in `work/` is ignored by git: nothing you do here is ever committed, and `git pull` never conflicts with it.

```
work/
    0_ADV_data_tank_v01.ipynb ... 3_SettlingFit_v01.ipynb   your copies of the notebooks
    data_OBS/       OBS records and step 1 outputs
    data_ADV/       ADV records: raw/, filtered/, figures/
    data_PSD/       particle size distribution csv files
    calibrations/   calibration files (step 2 writes here, step 3 reads from here)
    results/        step 3 results table and fit_plots/
```

All paths in the notebooks (`data_OBS/`, `calibrations/`, `results/`, ...) are relative to this folder. The notebooks start with the sample data (`MicroSphere_75_90_SAMPLE.csv`, `ADV_02_60rpm_SAMPLE.dat`), so they run as they are.

## Working with your own data

* Put raw records in `data_OBS/` and `data_ADV/raw/`, PSD files in `data_PSD/`, and change the inputs in the notebooks here.
* Keep as many copies of a notebook as you like, for example one step 1 notebook per record (`cp 1_OBS_timeseries_v04.ipynb 1_OBS_run07.ipynb`). Each saved copy is a record of the inputs you chose and the plots you chose them from. Copies in subfolders (e.g. `2026-11/`) work too.

## When the code is updated

* `git pull` in the project folder updates `tools/`, and your notebooks here use it the next time they run.
* Updated template notebooks do **not** replace your copies here. Run `python tools/setup_work.py` again: it adds new files but never overwrites, so to get a new version of a notebook you already have, rename or move your copy first, then rerun it (or copy the template from the project folder by hand).

To start over, delete or rename `work/` and run `python tools/setup_work.py`.
