"""Export the numbers the paper cites into `paper/numbers/*.json`, which `paper/paper_plan.typ` reads.

    python scripts/export_paper_numbers.py

Every `[[table]]` block in `paper/sources.toml` names a source CSV, a row filter, the key columns
and the cited values, each with a format. Writes one JSON per block: the source path and SHA-256,
and the matching rows with every value already formatted as the paper prints it. The paper typesets
those strings and never a number typed by hand, so a rerun reaches the text by re-exporting, and
`scripts/check_paper_provenance.py` re-exports in memory and fails on any difference.

Formats: `3f` / `2f` (fixed decimals), `pct1` / `pct0` (a fraction as a signed percentage with one
or no decimal, no sign when it rounds to zero, no `%`). Standard library only.
"""

from __future__ import annotations

import csv
import hashlib
import json
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
    if fmt in ("3f", "2f"):
        return f"{value:.{int(fmt[0])}f}"
    if fmt in ("pct1", "pct0"):
        return _signed(100 * value, int(fmt[-1]))
    raise ValueError(f"unknown format {fmt!r}")


def _matches(row: dict[str, str], where: dict) -> bool:
    return all(row[column] in wanted if isinstance(wanted, list) else row[column] == wanted
               for column, wanted in where.items())


def export_table(spec: dict) -> dict:
    """The JSON content of one `[[table]]` block, built from its source CSV."""
    source = ROOT / spec["source"]
    with open(source, newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if _matches(row, spec.get("where", {}))]
    if not rows:
        raise ValueError(f"{spec['name']}: no row of {spec['source']} matches {spec.get('where')}")
    exported = [{**{key: row[key] for key in spec["keys"]},
                 **{name: format_value(row[column], fmt) for name, (column, fmt) in spec["values"].items()}}
                for row in rows]
    return {"source": spec["source"], "source_sha256": file_sha256(source), "rows": exported}


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
