"""The notebooks of the RASPA workshop run against this raspalib.

The workshop notebooks are the reference for how raspalib is used: everything is defined in Python with keywords,
the phases of a driver are called one by one, and the results are read from the systems afterwards. Every code cell
of every notebook is executed here, with the cycle counts of the drivers capped so that a notebook takes seconds.

The notebooks are not part of this repository. Point RASPA_WORKSHOP_DIR at a checkout; without it the tests are
skipped. They need the packages of the workshop (numpy, scipy, matplotlib, tqdm, ase, py3Dmol, ipywidgets).

    RASPA_WORKSHOP_DIR=~/Raspa-Workshop-2026 python -m unittest discover -s tests/python -p "*_test.py"

RASPA_WORKSHOP_CYCLES sets the cap on the production cycles (default 100); the other phases get half of it.

With so few cycles the statistics are poor, and an analysis cell may be refused by numpy or scipy (a fit over an
empty time window, for instance). Such a cell is counted as skipped, and so is a later cell that needs a name the
skipped cell would have defined. An exception raised anywhere else, by raspalib or by the cell itself, fails the
notebook.
"""

import ast
import contextlib
import io
import json
import os
import pathlib
import sys
import tempfile
import traceback
import unittest

os.environ.setdefault("RASPALIB_QUIET", "1")
os.environ.setdefault("TQDM_DISABLE", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

WORKSHOP_DIRECTORY = os.environ.get("RASPA_WORKSHOP_DIR")
PRODUCTION_CYCLES = int(os.environ.get("RASPA_WORKSHOP_CYCLES", "100"))

# an exception that leaves one of these packages is the package refusing the data, not a fault of the notebook
NUMERICAL_PACKAGES = ("numpy", "scipy", "matplotlib")

CYCLE_KEYWORDS = {
    "number_of_cycles": PRODUCTION_CYCLES,
    "number_of_production_cycles": PRODUCTION_CYCLES,
    "number_of_pre_initialization_cycles": PRODUCTION_CYCLES // 2,
    "number_of_initialization_cycles": PRODUCTION_CYCLES // 2,
    "number_of_equilibration_cycles": PRODUCTION_CYCLES // 2,
}


def notebooks() -> list[pathlib.Path]:
    if not WORKSHOP_DIRECTORY:
        return []
    return sorted(pathlib.Path(WORKSHOP_DIRECTORY).expanduser().glob("day_*/*.ipynb"))


def code_cells(notebook: pathlib.Path) -> list[str]:
    cells = json.loads(notebook.read_text(encoding="utf-8"))["cells"]
    sources = []
    for cell in cells:
        if cell["cell_type"] != "code":
            continue
        lines = "".join(cell["source"]).splitlines()
        # the magics of IPython are not Python
        sources.append("\n".join(line for line in lines if not line.lstrip().startswith(("%", "!"))))
    return sources


def capped(driver):
    """The driver, constructed with no more cycles than the cap; the keywords are those the notebook wrote."""

    def construct(*arguments, **keywords):
        for keyword, cap in CYCLE_KEYWORDS.items():
            if keyword in keywords:
                keywords[keyword] = min(keywords[keyword], cap)
        return driver(*arguments, **keywords)

    return construct


def raised_inside_a_numerical_package(error: BaseException) -> bool:
    frames = traceback.extract_tb(error.__traceback__)
    parts = pathlib.Path(frames[-1].filename).parts if frames else ()
    return any(package in parts for package in NUMERICAL_PACKAGES)


def names_defined_by(source: str) -> set[str]:
    defined = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
            defined.add(node.id)
        elif isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            defined.add(node.name)
    return defined


def run_notebook(notebook: pathlib.Path) -> tuple[int, list[str]]:
    """Executes the code cells in one namespace; returns the number of cells run and the cells that were skipped."""
    import matplotlib.pyplot
    import raspalib

    namespace = {"__name__": "__main__"}
    missing: set[str] = set()
    skipped: list[str] = []
    cells = code_cells(notebook)

    originals = {name: getattr(raspalib, name) for name in ("MonteCarlo", "MolecularDynamics")}
    previous_directory = os.getcwd()
    with tempfile.TemporaryDirectory() as directory:
        os.chdir(directory)
        for name, driver in originals.items():
            setattr(raspalib, name, capped(driver))
        try:
            for index, source in enumerate(cells):
                label = f"{notebook.name} cell {index}"
                try:
                    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                        exec(compile(source, label, "exec"), namespace)
                except NameError as error:
                    if error.name not in missing:
                        raise AssertionError(f"{label}: {type(error).__name__}: {error}") from error
                    missing |= names_defined_by(source) - namespace.keys()
                    skipped.append(f"{label}: needs '{error.name}' of a skipped cell")
                except Exception as error:
                    if not raised_inside_a_numerical_package(error):
                        raise AssertionError(f"{label}: {type(error).__name__}: {error}") from error
                    missing |= names_defined_by(source) - namespace.keys()
                    skipped.append(f"{label}: {type(error).__name__}: {str(error).splitlines()[0] if str(error) else ''}")
                finally:
                    matplotlib.pyplot.close("all")
        finally:
            for name, driver in originals.items():
                setattr(raspalib, name, driver)
            os.chdir(previous_directory)
    return len(cells), skipped


@unittest.skipUnless(notebooks(), "RASPA_WORKSHOP_DIR does not point at the workshop notebooks")
class WorkshopNotebooks(unittest.TestCase):
    pass


def add_test(notebook: pathlib.Path) -> None:
    def test(self):
        number_of_cells, skipped = run_notebook(notebook)
        for line in skipped:
            print(f"  skipped: {line}", file=sys.stderr)
        # a notebook of which most cells are skipped has not been tested
        self.assertLess(len(skipped), max(2, number_of_cells // 2), f"too many cells skipped:\n" + "\n".join(skipped))

    name = "test_" + notebook.parent.name + "_" + notebook.stem
    test.__name__ = name
    setattr(WorkshopNotebooks, name, test)


for _notebook in notebooks():
    add_test(_notebook)


if __name__ == "__main__":
    unittest.main()
