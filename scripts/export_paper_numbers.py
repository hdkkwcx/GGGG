"""Export the numbers the paper cites into `paper/numbers/*.json`, which `paper/paper_plan.typ` reads.

    python scripts/export_paper_numbers.py

Every `[[table]]` block in `paper/sources.toml` names a source CSV, a row filter, the key columns
and the cited values, each with a format. Writes one JSON per block: the source path and SHA-256,
and the matching rows with every value already formatted as the paper prints it. A block's optional
`[table.stats.<name>]` entries add the numbers worked out from several rows (a range, a median, a
largest difference) under `stats`, so the paper looks those up too and computes nothing. The paper typesets
those strings and never a number typed by hand, so a rerun reaches the text by re-exporting, and
`scripts/check_paper_provenance.py` re-exports in memory and fails on any difference.

Formats: `Nf` (N fixed decimals, e.g. `3f`, `0f`), `+Nf` (the same with a sign, for a value the text
prints signed, such as a bias or a value already in percent), `pct1` / `pct0` (a fraction as a signed
percentage with one or no decimal, no `%`). A signed format drops the sign when the value rounds to
zero. Standard library only.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ROOT / "paper" / "sources.toml"
DATA_DIR = ROOT / "paper" / "numbers"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _signed(value: float, decimals: int) -> str:
    text = f"{abs(value):.{decimals}f}"
    if float(text) == 0:
        return text
    return ("+" if value > 0 else "-") + text


def format_value(raw: str, fmt: str) -> str:
    value = float(raw)
    if fmt in ("pct1", "pct0"):
        return _signed(100 * value, int(fmt[-1]))
    if fmt.endswith("f") and fmt.lstrip("+")[:-1].isdigit():
        decimals = int(fmt.lstrip("+")[:-1])
        return _signed(value, decimals) if fmt.startswith("+") else f"{value:.{decimals}f}"
    raise ValueError(f"unknown format {fmt!r}")


def _matches(row: dict[str, str], where: dict) -> bool:
    return all(row[column] in wanted if isinstance(wanted, list) else row[column] == wanted
               for column, wanted in where.items())


def compute_stat(rows: list[dict[str, str]], name: str, spec: dict) -> str:
    """One `[table.stats.<name>]` entry: a statistic over the matching rows of one column, formatted.

    `stat` is `min`, `max`, `median`, `range` ("min–max") or `maxgap`, the largest |a − b| between
    levels `a` and `b` of column `key` over rows that agree on the `by` columns. `abs = true` takes
    the size of each value first. Computed from the unrounded values; a leading `+` is dropped.
    """
    keep = [row for row in rows if _matches(row, spec.get("where", {}))]
    column, fmt = spec["column"], spec["fmt"]

    def show(value: float) -> str:
        return format_value(repr(value), fmt).lstrip("+")

    if spec["stat"] == "maxgap":
        partner = {tuple(row[k] for k in spec["by"]): row for row in keep if row[spec["key"]] == spec["b"]}
        values = [abs(float(row[column]) - float(partner[pair][column]))
                  for row in keep if row[spec["key"]] == spec["a"]
                  and (pair := tuple(row[k] for k in spec["by"])) in partner]
        spec = {**spec, "stat": "max"}
    else:
        values = [float(row[column]) for row in keep]
        values = [abs(value) if spec.get("abs") else value for value in values if not math.isnan(value)]
    if not values:
        raise ValueError(f"stat {name}: no value of {column} matches {spec.get('where')}")
    if spec["stat"] == "range":
        return f"{show(min(values))}–{show(max(values))}"
    return show({"min": min, "max": max, "median": statistics.median}[spec["stat"]](values))


def export_table(spec: dict) -> dict:
    """The JSON content of one `[[table]]` block, built from its source CSV."""
    source = ROOT / spec["source"]
    with open(source, newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if _matches(row, spec.get("where", {}))]
    if not rows:
        raise ValueError(f"{spec['name']}: no row of {spec['source']} matches {spec.get('where')}")
    exported = [{**{key: row[key] for key in spec["keys"]},
                 **{name: format_value(row[column], fmt) for name, (column, fmt) in spec["values"].items()}}
                for row in rows] if "values" in spec else []
    content = {"source": spec["source"], "source_sha256": file_sha256(source), "rows": exported}
    if "stats" in spec:
        content["stats"] = {name: compute_stat(rows, name, stat) for name, stat in spec["stats"].items()}
    return content


def table_specs() -> list[dict]:
    return tomllib.loads(SOURCES.read_text(encoding="utf-8")).get("table", [])


def render(content: dict) -> str:
    return json.dumps(content, indent=1, ensure_ascii=False) + "\n"


def main() -> int:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for spec in table_specs():
        path = DATA_DIR / f"{spec['name']}.json"
        path.write_text(render(export_table(spec)), encoding="utf-8")
        print(f"wrote {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
