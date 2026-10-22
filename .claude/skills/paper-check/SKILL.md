---
name: paper-check
description: Check, without changing anything, that the paper is consistent and up to date - every figure and cited number traced to committed code, every make_paper.jl step reusable, and paper_plan.typ compiling cleanly. Use when the user asks to check the paper, before committing changes under paper/, or after analysis code changed.
---

# Paper check

Read-only. Run from the repository root and report; never edit, commit or rerun anything here.

1. **Provenance.** `"${MIXEDGWR_PYTHON:-python3.13}" scripts/check_paper_provenance.py`
   It must end with `Every entry traces to committed code.` Report every `FAIL` line as it is.
   (`python3` on this Mac is 3.9 and lacks `tomllib`; use `MIXEDGWR_PYTHON` or `python3.13`.)

2. **Steps.** `julia --project=. scripts/make_paper.jl` (the default is a dry run).
   Every analysis and figure step should be `REUSE`; data steps `OK` or `ABSENT`. For any other:
   - `CHANGED`: quote the files it lists. The fix is the user's call: review the diff
     (`git diff <commit> -- <files>`), then either rerun the step (`--run --force <step>`) or bump
     the step's `reviewed` commit in `paper/sources.toml` if the change cannot move the result.
   - `STALE` / `MISSING`: quote the reason and the command the dry run prints.
   Never pass `--run` yourself: it reruns analyses (the benchmark alone takes about 4.3 h) and
   needs the user's go-ahead.

3. **Compile.** `typst compile paper/paper_plan.typ` must finish with no warning or error.
   To read the text, `pdftotext paper/paper_plan.pdf -`.

Report a short table: provenance (OK or the failures), steps (counts by status and any non-OK
step), compile (clean or the messages). If everything passes, say so in one line.
