# CLAUDE.md

Guidance for Claude Code sessions in this repo. The README covers what the toolbox does and how users run it; this file covers how to work on it.

## Layout in one paragraph

Template notebooks (`0_…` to `3_…ipynb`) sit at the top level, shared code is in `tools/`, and docs are in `docs/` (`01_theory.md`, `02_user_guide.md`, `03_reference.md`). Users never work in the templates: `python tools/setup_work.py` copies them plus `work_template/` (sample records, example PSDs, example calibration) into `work/`, which is git-ignored and holds all data, calibrations and results. Every data path in a notebook is relative to `work/`.

## Rules

- **Never commit anything from `work/`**, and don't add data or outputs elsewhere in the tree. Only the files in `work_template/` are tracked sample data.
- **Edit the templates, not the `work/` copies,** when changing a notebook for everyone. `setup_work.py` never overwrites, so changed templates reach `work/` only under a new name; say so when a template changes.
- **Keep the first code cell of each notebook identical** across notebooks: it finds the project folder (the one holding `tools/` and `work_template/`), puts `tools/` on `sys.path`, runs `setup_work` if `work/` is missing, and `chdir`s into `work/`.
- **Put processing logic in `tools/`, not in notebook cells.** `git pull` updates `tools/` for everyone's `work/` copies, but never their notebooks. Notebooks should hold inputs, calls and plots.
- **Notebooks are committed without outputs.** nbstripout runs as a git filter (`nbstripout --install` once per clone). Don't commit a notebook in a clone where it isn't installed.
- **Keep docs in step with code.** Changes to inputs, file names, outputs or folders usually touch `README.md` and all three files in `docs/`; check them before finishing.

## Naming conventions

- The step 1 output is `*_time_series_decayONLY.csv`, held in the variable `OBS_decay` (formerly `OBS_after_paddle`; don't reintroduce it).
- Say "mixer", not "paddle", in text and comments.
- Sample inputs carry `_SAMPLE` in their names (`MicroSphere_75_90_SAMPLE.csv`, `ADV_02_60rpm_SAMPLE.dat`). Notebooks 0 and 1 pick their default file by name with `filenames.index(...)`, not by position.
- Folder names are lower case: `data_OBS/`, `data_ADV/{raw,filtered,figures}/`, `data_PSD/`, `calibrations/`, `results/`.

## Testing a change

There is no test suite. To check a change end to end, without touching the real `work/`:

1. Copy only what git tracks into a scratch folder (for example `git archive HEAD | tar -x -C <scratch>`, plus any uncommitted edits).
2. In that copy run `python tools/setup_work.py`, then execute the notebooks from `work/` in order 1 → 2 → 3 (and 0 if ADV code changed) with `jupyter nbconvert --to notebook --execute`.
3. Expected on the sample: step 2 reproduces `work_template/calibrations/cal_60rpm.json` (K₀, K_min, rms), and step 3 gives a density of about 1199 kg m⁻³ (known: 1200).

Notebooks use `%matplotlib widget` (ipympl), so steps 0 and 1 are interactive; values such as the averaging windows and `mixer_off` are read off those plots by the user.

## Planned direction

The `work/` copy approach is a first step. The intended next move is per-notebook **settings files** (inputs in git-ignored files kept in `work/`, with tracked examples that point at the SAMPLE data), with the inline processing in notebooks 0 and 1 moved into `tools/` first. Propose this when touching how notebooks take their inputs.
