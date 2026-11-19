"""
setup_work: create (or top up) the work/ folder in which the notebooks are run.

    python tools/setup_work.py        # run once from the project folder after cloning

Copies into work/:
  * everything in work_template/ (sample OBS and ADV records, example PSD files, example calibration, README),
  * the template notebooks (0_*.ipynb ... 3_*.ipynb) from the project folder,
and creates work/results/.

Existing files in work/ are NEVER overwritten, so it is safe to run again at any time, for example after
`git pull` to pick up new sample files or new versions of the notebooks (a template whose name already exists
in work/ is left alone; copy it by hand under a new name if you want the updated version).

The notebooks call setup_work() themselves if work/ does not exist yet.
"""
import glob
import os
import shutil

SKIP = {".DS_Store"}


def project_root():
    """The project folder: the parent of the tools/ folder this file lives in."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def setup_work(root=None, verbose=True):
    """Create or top up <root>/work from <root>/work_template and the template notebooks.
    Returns the list of files added (paths relative to the project folder)."""
    root = root or project_root()
    src = os.path.join(root, "work_template")
    dst = os.path.join(root, "work")
    if not os.path.isdir(src):
        raise FileNotFoundError(f"{src} not found; is {root} the project folder?")

    added = []

    def copy_new(source, target):
        if os.path.exists(target):
            return
        os.makedirs(os.path.dirname(target), exist_ok=True)
        shutil.copy2(source, target)
        added.append(os.path.relpath(target, root))

    # 1. the starter files (samples, example calibration, README)
    for dirpath, _, filenames in os.walk(src):
        rel = os.path.relpath(dirpath, src)
        for name in filenames:
            if name not in SKIP:
                copy_new(os.path.join(dirpath, name), os.path.normpath(os.path.join(dst, rel, name)))

    # 2. the template notebooks
    for nb in sorted(glob.glob(os.path.join(root, "[0-9]_*.ipynb"))):
        copy_new(nb, os.path.join(dst, os.path.basename(nb)))

    # 3. the output folder of step 3
    os.makedirs(os.path.join(dst, "results"), exist_ok=True)

    if verbose:
        if added:
            print(f"work/ set up in {dst}; added {len(added)} file(s):")
            for a in added:
                print("   ", a)
        else:
            print(f"work/ is up to date ({dst}); nothing added")
    return added


if __name__ == "__main__":
    setup_work()
