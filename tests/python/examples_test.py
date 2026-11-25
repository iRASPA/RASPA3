"""The run.py scripts of the examples run against this raspalib.

Every example directory has the input of the command-line program (simulation.json with the force field, the
molecules and the framework next to it) and a run.py that sets up the same simulation with raspalib, reading the
force field and the molecules from those files. Each script is executed here in a copy of its directory, with the
cycle counts of the drivers capped so that it takes seconds.

    python -m unittest discover -s tests/python -p "*_test.py"

RASPA_EXAMPLES_CYCLES sets the cap on the production cycles (default 20); the other phases get half of it.
"""

import contextlib
import io
import os
import pathlib
import runpy
import shutil
import subprocess
import tempfile
import unittest

os.environ.setdefault("RASPALIB_QUIET", "1")

import raspalib

ROOT = pathlib.Path(__file__).resolve().parents[2]
PRODUCTION_CYCLES = int(os.environ.get("RASPA_EXAMPLES_CYCLES", "20"))

CYCLE_KEYWORDS = {
    "number_of_cycles": PRODUCTION_CYCLES,
    "number_of_production_cycles": PRODUCTION_CYCLES,
    "number_of_pre_initialization_cycles": PRODUCTION_CYCLES // 2,
    "number_of_initialization_cycles": PRODUCTION_CYCLES // 2,
    "number_of_equilibration_cycles": PRODUCTION_CYCLES // 2,
}


def scripts() -> list[pathlib.Path]:
    """The scripts of the repository that use raspalib (those under version control, when git can tell)."""
    try:
        listed = subprocess.run(
            ["git", "ls-files", "examples/**/run.py"], cwd=ROOT, capture_output=True, text=True, check=True
        ).stdout.split()
        found = [ROOT / name for name in listed]
    except (OSError, subprocess.CalledProcessError):
        found = list((ROOT / "examples").glob("**/run.py"))
    return sorted(script for script in found if "import raspalib" in script.read_text(encoding="utf-8"))


def capped(driver):
    def construct(*arguments, **keywords):
        for keyword, cap in CYCLE_KEYWORDS.items():
            if keyword in keywords:
                keywords[keyword] = min(keywords[keyword], cap)
        return driver(*arguments, **keywords)

    return construct


def run_script(script: pathlib.Path) -> str:
    originals = {name: getattr(raspalib, name) for name in ("MonteCarlo", "MolecularDynamics")}
    previous_directory = os.getcwd()
    output = io.StringIO()
    with tempfile.TemporaryDirectory() as directory:
        # the input files only: an example directory may hold the output of an earlier run
        for item in script.parent.iterdir():
            if item.is_file() and item.suffix in (".json", ".cif", ".py"):
                shutil.copy(item, directory)
        os.chdir(directory)
        for name, driver in originals.items():
            setattr(raspalib, name, capped(driver))
        try:
            with contextlib.redirect_stdout(output):
                runpy.run_path(str(pathlib.Path(directory) / "run.py"), run_name="__main__")
        finally:
            for name, driver in originals.items():
                setattr(raspalib, name, driver)
            os.chdir(previous_directory)
    return output.getvalue()


class Examples(unittest.TestCase):
    def test_there_are_examples(self):
        self.assertGreater(len(scripts()), 20)


def add_test(script: pathlib.Path) -> None:
    def test(self):
        output = run_script(script)
        self.assertIn("average energy", output)

    name = "test_" + script.parent.parent.name + "_" + script.parent.name
    test.__name__ = name
    setattr(Examples, name, test)


for _script in scripts():
    add_test(_script)


if __name__ == "__main__":
    unittest.main()
