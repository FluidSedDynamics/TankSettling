# Plan: move notebook inputs into settings files

Branch: `settings-files` · Status: **draft for review, no code changed yet** · Delete this file before merging into `main`.

## Goal

Today every user works in copies of the notebooks inside `work/`. That keeps tracked files clean, but `git pull` never reaches those copies: a fixed or improved notebook only arrives as a new file, and the user has to move their inputs into it by hand.

After this change:

- Everyone runs the **tracked notebooks** in the project folder. They hold code only and are updated by `git pull`.
- Each user's inputs live in **settings files** in `work/settings/`, which git ignores. Data, calibrations and results stay in `work/` as now.
- Running a notebook never changes a tracked file in a way git sees (nbstripout removes outputs; see §6).

## 1. Settings format: TOML

One file per step, in `work/settings/`:

```
work/settings/
    step0_adv.toml
    step1_obs.toml
    step2_calibration.toml
    step3_fit.toml
```

Why TOML: it allows comments (so every setting can carry its explanation, as the input cells do now), it is easy to read and edit by hand, and Python reads it with the standard library (`tomllib`, Python 3.11+), so no new dependency.

TOML has no `None`. Convention: **leave a setting out (or comment it out) to get the automatic/default value.** Every default lives in code in `tools/`, so a settings file only needs the values a user actually changes, and a newer version of the code can add settings without breaking older files.

## 2. One file per step, with named entries

Steps 0 and 1 are run once per record, and the values chosen for each record are worth keeping. So each file holds a **selector** at the top and one **named entry** per record (or calibration, or batch of runs):

```toml
# work/settings/step1_obs.toml
record = "MicroSphere_75_90_SAMPLE"     # which record the notebook processes now

[defaults]                              # apply to every record unless overridden below
time_correction = { hours = 0, minutes = -1, seconds = 2 }
move_third_column_to_end = true

[records.MicroSphere_75_90_SAMPLE]
Co_window_s         = [120, 340]        # averaging window for the initial concentration
background_window_s = [700, 1014]       # averaging window for the background
mixer_off_s         = 350               # time the mixer was switched off

# optional NTU -> mg/L conversion (single-size particle standards)
mg_per_L = { slope = 54.1703, intercept = -87.3999 }
```

The same pattern for the other steps:

| File | Selector | Named entries | Contents |
|---|---|---|---|
| `step0_adv.toml` | `record` | `[records.<file name>]` | `fs`, despike options, `mixer_off_s`, `last_sec_for_avg` |
| `step1_obs.toml` | `record` | `[records.<file name>]` | time correction, averaging windows, mixer off, optional mg/L |
| `step2_calibration.toml` | `calibration` | `[calibrations.<name>]` | files, rpm, geometry, water, standard, turbulence model, ADV shape, validation |
| `step3_fit.toml` | `batch` | `[batches.<name>]` with `[[batches.<name>.runs]]` | calibration file, defaults, the list of runs, output names |

The step 3 run list maps directly onto TOML's array of tables:

```toml
batch = "sample"

[batches.sample]
calibration_file = "calibrations/cal_60rpm.json"
rpm = 60
what_to_fit = "rhos"

[[batches.sample.runs]]
name     = "standard 75-90 um (PSD file)"
obs_file = "data_OBS/MicroSphere_75_90_SAMPLE_time_series_decayONLY.csv"
psd_file = "data_PSD/PSD_example_75_90.csv"

[[batches.sample.runs]]
name     = "standard, single size 82.5 um"
obs_file = "data_OBS/MicroSphere_75_90_SAMPLE_time_series_decayONLY.csv"
psd_file = "data_PSD/PSD_example_single_82_5um.csv"
```

## 3. The interactive steps (0 and 1)

Steps 0 and 1 need values read off plots. The loop becomes: run the notebook, read the windows off the plot, type them into `step1_obs.toml` (open beside the notebook in VS Code or Jupyter Lab), rerun the cell. The notebook prints the path of the settings file it used and the values it read, so it is always clear what was applied.

A new record with no entry yet still runs as far as the first plot (needed to choose the windows), then stops with a message giving the entry to add.

Later, optional: pick the windows by clicking on the ipympl plot and have the notebook write them into the settings file. Not part of this branch.

## 4. Keeping a record of each run

With working copies, a saved notebook recorded the values *and* the plots for each record. Notebooks will now be rerun for each record, so that record has to come from the outputs instead:

- **Step 0 and step 1:** save the plots to `work/data_ADV/figures/` and `work/data_OBS/figures/`, and write the settings used next to each output (`<record>_settings_used.toml`).
- **Step 2:** the calibration JSON already stores its inputs; add any that are missing.
- **Step 3:** write the batch's settings next to the results csv (`results/<batch>_settings_used.toml`).
- **Optional:** a setting to save an executed copy of the notebook (HTML) into `work/` for runs worth archiving.

## 5. Code changes

**New `tools/settings.py`:**
- `load_settings(step)`: finds `work/settings/<file>`, copies it from `work_template/settings/` if missing, reads it, merges `[defaults]` and the selected entry over the code defaults, and returns them.
- **Rejects unknown setting names** with a message suggesting the closest valid name, so a typo can't be silently ignored.
- `save_settings_used(...)`: writes the record described in §4.

**Processing moved out of notebooks 0 and 1 into `tools/`** (new `tools/obs_tools.py` and `tools/adv_tools.py`): reading the OBS logger file and reshaping it into a time series, window averages, extracting the decay; reading Vectrino `.dat` files, despiking, extracting the decay. Notebooks 2 and 3 already call `settling_tools` for nearly everything and only need their input cells replaced.

**Notebooks:** each input cell becomes `cfg = load_settings('step1')` (or the matching step). The first cell (find the project folder, `chdir` into `work/`) stays as it is.

**`tools/setup_work.py`:** creates `work/settings/` from `work_template/settings/` and **stops copying notebooks** into `work/`. It also checks that nbstripout is installed for the clone and installs it if not (see §6).

**`work_template/settings/`:** the four example settings files, set up for the SAMPLE data, so a fresh clone runs end to end with no edits.

Changes to `settling_tools.py` stay additive, so existing `work/` notebook copies (such as the real ones on `main`) keep running against the new code.

## 6. nbstripout becomes required

Running a tracked notebook writes outputs into it. With nbstripout installed, git ignores the outputs and `git status` stays clean. Without it, every run shows the notebook as modified, and `git pull` can refuse to update it.

So: add `nbstripout` to `requirements.txt`, have `setup_work.py` run `nbstripout --install` for the clone, and have the first cell of each notebook warns if the filter isn't configured.

## 7. Order of work and checks

Each phase ends with the sample run end to end in this worktree (`work/` here holds SAMPLE data only) and outputs compared with `main`:

1. **Move notebook 0 and 1 processing into `tools/`**, notebooks still using inline inputs. Check: `*_filtered_decay_only.csv` and `*_time_series_decayONLY.csv` identical to `main`.
2. **Add `tools/settings.py`** and the four example files. Check: loading each example gives exactly the values in today's input cells.
3. **Convert the notebooks one at a time** (1, 0, 2, 3). Check after each: same outputs as `main`; step 2 reproduces `cal_60rpm.json` (K₀, K_min, rms); step 3 gives about 1199 kg m⁻³.
4. **Update `setup_work.py`** and test on a fresh copy of the tracked files: setup, run all four notebooks in place, confirm `git status` is clean.
5. **Docs:** README (quick start, workflow, layout), `docs/02_user_guide.md` (installation, every step's inputs), `docs/03_reference.md` (every setting, with defaults), `work_template/README.md`, `CLAUDE.md`.
6. **Pull request** into `main`, then migrate the real `work/` on `main` (§8). Delete this plan file before merging.

## 8. Moving existing work over

For anyone with a `work/` from the current layout (including the real one on `main`):

- Nothing in `work/` is deleted or overwritten; old notebook copies keep working.
- `python tools/setup_work.py` adds `work/settings/` with the example files.
- Values from old notebook copies are copied into the settings files by hand, or with a small helper script that reads a `work/` notebook's input cells and writes the matching settings entry. That helper is worth it only if there are many records to move; decide in phase 6.

## Open questions

1. **Named entries in one file per step (proposed), or one settings file per record?** One file per step keeps everything for a step in one place and needs no new files per record; one file per record is simpler for each file but adds many small files.
2. **Run the tracked notebooks in place (proposed), or keep copying notebooks into `work/` and only move the inputs out?** In place means `git pull` updates notebooks for everyone; copying keeps the current habit but leaves the update problem unsolved.
3. **Where to record what was run (§4):** is a `_settings_used.toml` next to each output enough, or should saving an executed HTML copy be on by default?
4. **Step 1's optional mg/L conversion:** keep it in step 1's settings, or drop it from the workflow since the fit doesn't use it?
