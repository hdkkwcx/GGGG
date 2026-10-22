"""Check that every figure and table in the paper traces back to committed code.

    python scripts/check_paper_provenance.py

Reads `paper/sources.toml` and, for every entry, checks that:

- the source output exists, and the paper's copy of a figure is byte-identical to it;
- every provenance file it lists records `git_dirty=false` (and `git_dirty_at_end` /
  `report_git_dirty=false` where present), so the commit it names reproduces it;
- every `input_sha256:<path>` row a run recorded still matches that file under `data/`. A run from
  before inputs were hashed gets a note, not a failure: its code is still pinned.

Every `[[table]]` block, the numbers the text cites, is also re-exported in memory
(`export_paper_numbers.py`) and must equal its tracked `paper/numbers/<name>.json`, so a hand edit,
a stale export or a changed source all fail.

It also fails when `paper/figures/` holds a file, or `paper/paper_plan.typ` includes an image, that
has no entry, or when an entry names no `step`, or one that is not a `[[step]]` block (the graph
`scripts/make_paper.jl` rebuilds the paper from). One line per entry, then the commits involved;
the exit status is nonzero if anything failed. Standard library only.
"""

from __future__ import annotations

import csv
import hashlib
import re
import sys
import tomllib
from pathlib import Path

from export_paper_numbers import export_table, render

ROOT = Path(__file__).resolve().parents[1]
CLEAN_KEYS = ("git_dirty", "git_dirty_at_end", "report_git_dirty")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_provenance(path: Path) -> dict[str, str]:
    with open(path, newline="", encoding="utf-8") as handle:
        return {row["key"]: row["value"] for row in csv.DictReader(handle)}


def check_item(item: dict, commits: dict[str, set[str]]) -> tuple[list[str], list[str]]:
    """The entry's failures, and its notes: gaps that do not make it untraceable."""
    problems, notes = [], []
    source = ROOT / item["source"]
    if not source.exists():
        return [f"source missing: {item['source']}"], notes
    if "paper" in item:
        paper = ROOT / item["paper"]
        if not paper.is_file():
            problems.append(f"paper file missing: {item['paper']}")
        elif sha256(paper) != sha256(source):
            problems.append("paper copy differs from its source")
    for name in item["provenance"]:
        path = ROOT / name
        if not path.is_file():
            problems.append(f"no provenance: {name}")
            continue
        record = read_provenance(path)
        short = path.parent.name + "/" + path.name
        for key in CLEAN_KEYS:
            if key in record and record[key] != "false":
                problems.append(f"{key}={record[key]} in {short}")
        for key in ("git_commit", "report_git_commit"):
            if key in record:
                commits.setdefault(record[key][:7], set()).add(short)
        inputs = {key.removeprefix("input_sha256:"): value for key, value in record.items()
                  if key.startswith("input_sha256:")}
        for relative, recorded in inputs.items():
            data = Path(relative) if Path(relative).is_absolute() else ROOT / relative
            if not data.is_file():
                problems.append(f"input gone: {relative}")
            elif sha256(data) != recorded:
                problems.append(f"input changed since the run: {relative}")
        if name.endswith(("benchmark_scope.csv", "/provenance.csv")) and not inputs:
            notes.append(f"input checksums not recorded in {short}")
    return problems, notes


def unlisted(items: list[dict]) -> list[str]:
    listed = {item["paper"] for item in items if "paper" in item}
    on_disk = {f"paper/figures/{path.name}" for path in (ROOT / "paper" / "figures").glob("*")
               if path.is_file() and not path.name.startswith(".")}
    included = {f"paper/{name}" for name in
                re.findall(r'image\("([^"]+)"', (ROOT / "paper" / "paper_plan.typ").read_text(encoding="utf-8"))}
    return sorted((on_disk | included) - listed)


def check_table(spec: dict, commits: dict[str, set[str]]) -> tuple[list[str], list[str]]:
    """A `[[table]]` block: the source's provenance, then the export against the tracked JSON."""
    problems, notes = check_item(spec, commits)
    if problems:
        return problems, notes
    tracked = ROOT / "paper" / "numbers" / f"{spec['name']}.json"
    if not tracked.is_file():
        problems.append(f"not exported: paper/numbers/{spec['name']}.json")
    elif tracked.read_text(encoding="utf-8") != render(export_table(spec)):
        problems.append(f"paper/numbers/{spec['name']}.json differs from a fresh export of its source")
    return problems, notes


def main() -> int:
    sources = tomllib.loads((ROOT / "paper" / "sources.toml").read_text(encoding="utf-8"))
    items, tables = sources["item"], sources.get("table", [])
    steps = {step["id"] for step in sources.get("step", [])}
    commits: dict[str, set[str]] = {}
    failed = 0
    width = max(len(entry["label"]) for entry in items + tables)
    for item in items + tables:
        problems, notes = (check_table if item in tables else check_item)(item, commits)
        if item.get("step") not in steps:
            problems.append(f"no step: {item['step']}" if "step" in item else "names no step")
        failed += bool(problems)
        status = "FAIL: " + "; ".join(problems) if problems else "OK"
        print(f"{item['label']:<{width}}  {status}" + "".join(f"\n{'':<{width}}  note: {note}" for note in notes))
    for path in unlisted(items):
        failed += 1
        print(f"{path:<{width}}  FAIL: no entry in paper/sources.toml")
    print("\nCommits: " + (", ".join(f"{commit} ({len(files)} files)" for commit, files in sorted(commits.items()))
                           or "none recorded"))
    print(f"{failed} of {len(items) + len(tables)} entries failed" if failed else "Every entry traces to committed code.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
