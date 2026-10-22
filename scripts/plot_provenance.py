"""Which code and libraries drew a figure directory: `plot_provenance.csv`, written beside the figures.

The Julia run's own provenance (`provenance.csv`, `benchmark_scope.csv`) describes the analysis, not
the plotting script, and the two can be different commits. Imported by the paper's plot scripts,
which run as `python scripts/<name>.py` and so have `scripts/` on `sys.path`.

Rows: the git commit and dirty flag (read when the figures are written), the command that drew them
(the script, relative to the repository, and its arguments), the Python, matplotlib and numpy
versions, and the SHA-256 of `scripts/requirements-plot.txt`, the environment they are pinned
to. A figure drawn straight from `data/` rather than from a Julia run's outputs passes those files as
`inputs`, which adds an `input_sha256:<path>` row for each, as `start_provenance` writes them. A dirty
tree is recorded and warned about, never refused, as in `src/util/RunProvenance.jl`.
"""

from __future__ import annotations

import csv
import hashlib
import platform
import subprocess
import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np


def _git(root: Path, *args: str) -> str:
    try:
        return subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def _command(root: Path) -> str:
    """The running script, relative to `root` when it lies under it, then its arguments."""
    script = Path(sys.argv[0]).resolve()
    try:
        script = script.relative_to(root.resolve())
    except ValueError:
        pass
    return " ".join([script.as_posix(), *sys.argv[1:]])


def write_plot_provenance(fig_dir: Path, root: Path, inputs: tuple[Path, ...] = ()) -> None:
    porcelain = _git(root, "status", "--porcelain")
    dirty = porcelain if porcelain == "unavailable" else str(porcelain != "").lower()
    if dirty == "true":
        print("warning: the working tree has uncommitted changes; plot_provenance.csv records "
              "git_dirty=true, and no commit reproduces these figures", file=sys.stderr)
    requirements = root / "scripts" / "requirements-plot.txt"
    rows = [["git_commit", _git(root, "rev-parse", "HEAD")], ["git_dirty", dirty], ["command", _command(root)],
            ["python", platform.python_version()], ["matplotlib", mpl.__version__], ["numpy", np.__version__],
            ["requirements_sha256",
             hashlib.sha256(requirements.read_bytes()).hexdigest() if requirements.is_file() else "missing"]]
    rows += [[f"input_sha256:{path.relative_to(root).as_posix()}", hashlib.sha256(path.read_bytes()).hexdigest()]
             for path in inputs]
    fig_dir.mkdir(parents=True, exist_ok=True)
    with open(fig_dir / "plot_provenance.csv", "w", newline="", encoding="utf-8") as handle:
        csv.writer(handle).writerows([["key", "value"], *rows])
