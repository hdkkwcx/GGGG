#!/usr/bin/env julia

# Do the paper's verdicts survive resampling whole weeks instead of days? Most of its 95% intervals
# and p-values come from a day-block bootstrap, which treats days as independent. FY4B's worst errors
# come in multi-day episodes (`fy4b_quality_handover.md`), so this recomputes them with 7- and 14-day
# blocks, from the same stored predictions, and lists every verdict that changes. A check: no paper
# text or number changes.
#
#   julia --project=. scripts/compare_interval_blocks.jl
#
# It reads the published (day-block) files and the same files written with `--block-days N` (the
# `_block{N}d` folders) by these, in order:
#     julia --project=. scripts/run_claim_reassessment.jl <paper run> --block-days N
#     julia --project=. scripts/run_gauge_satellite_fusion.jl <paper run> --block-days N
#     julia --project=. scripts/compare_fy4b_quality.jl --block-days N
# for N = 1 (the guard), 7 and 14.
#
# Families, and what makes a verdict in each:
# - RQ3 `fusion_vs_gauge_only.csv` (the `vs_adw_uncal` table and the text's "significant"): Holm p
#   across the stratum's family < 0.05. The figure's filled cells also need the relative interval to
#   exclude 0, so a row flips if either changes. Rows the paper cites: methods mgwr and the agreement
#   blend.
# - RQ3 `interpolation_ranking.csv`, `anchor_comparison.csv`, `detection_skill_bootstrap.csv`
#   (figures and the report's tables): the interval excludes 0.
# - §3.4 `fy4b_quality/paired.csv`: the RMSE interval excludes 0 (the text's `dq`, and "spans zero").
#   Rows the paper cites: those §3.4 prints.
#
# The reading, fixed before the runs:
# - A verdict *flips* when it differs between day blocks and 7-day blocks.
# - 14-day blocks are reported as a sensitivity only.
# - If no cited verdict flips at 7 days, the paper's verdicts stand.
# - If any does, the flipped rows are listed, and the user decides the wording.
#
# Checks; the script errors if one fails:
# - every `_block1d` file is byte-identical to its published file;
# - every point estimate is the same across block lengths, to a relative 1e-9 (summation order), so
#   only intervals and p-values move.
#
# Writes output/interval_block_check/: verdicts.csv (every row, with its intervals, p-values and
# verdicts per block length, whether the paper cites it, and how much its interval widens),
# summary.csv (per family) and provenance.csv.

const ROOT = normpath(joinpath(@__DIR__, ".."))

using CSV, DataFrames, Printf, Statistics

include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("RunProvenance", "TableIO")
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic

const RUN = "interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only_satwetblend_" *
    "blendagrenv_fusedanchor_nofy4b_gpmuncal"
const CLAIMS_DIR = joinpath(ROOT, "output", "benchmark_diagnostics", RUN)
const FUSION_DIR = joinpath(ROOT, "output", "gauge_satellite_fusion_nofy4b_gpmuncal")
const QUALITY_DIR = joinpath(ROOT, "output", "fy4b_quality")
const OUTDIR = joinpath(ROOT, "output", "interval_block_check")
const BLOCKS = [7, 14]
const JUDGED = 7
const ALPHA = 0.05
const CITED_RQ3_METHODS = ("mgwr", "blend_agrenv_mgwr")
"""§3.4's `dq` calls: (scheme, baseline, method)."""
const CITED_QUALITY = Set([
    ("balanced_spatial", "FY4B", "raw"), ("balanced_spatial", "FY4B", "mgwr"),
    ("balanced_spatial", "FY4B", "blend_agrenv_mgwr"), ("balanced_spatial", "GPM_uncal", "mgwr"),
    ("balanced_spatial", "GPM_uncal", "blend_agrenv_mgwr"), ("balanced_spatial", "GSMaP", "mgwr"),
    ("balanced_spatial", "GSMaP", "blend_agrenv_mgwr"), ("random", "FY4B", "mgwr"),
])

block_dir(dir, days) = dir * "_block$(days)d"
excludes_zero(lo, hi) = isfinite(lo) && isfinite(hi) && (lo > 0 || hi < 0)

"""
One family: which file, the key columns, the point/interval/p columns, how a row's verdict is made,
and whether the paper cites a row.
"""
const FAMILIES = [
    (; name="rq3_fusion_vs_adw", dir=FUSION_DIR, file="fusion_vs_gauge_only.csv",
        keys=[:product, :method, :stratum], point=:relative_improvement, lo=:relative_ci_low,
        hi=:relative_ci_high, p=:pvalue_holm_family,
        verdict=row -> (row.p < ALPHA, row.p < ALPHA && excludes_zero(row.lo, row.hi)),
        cited=row -> row.method in CITED_RQ3_METHODS),
    (; name="rq3_interpolation_ranking", dir=FUSION_DIR, file="interpolation_ranking.csv",
        keys=[:method, :stratum, :baseline], point=:relative_improvement, lo=:ci_low, hi=:ci_high,
        p=:pvalue_holm, verdict=row -> (excludes_zero(row.lo, row.hi),), cited=row -> false),
    (; name="rq3_anchor_comparison", dir=FUSION_DIR, file="anchor_comparison.csv",
        keys=[:method, :product], point=:relative_improvement, lo=:relative_ci_low,
        hi=:relative_ci_high, p=:pvalue_holm, verdict=row -> (excludes_zero(row.lo, row.hi),),
        cited=row -> false),
    (; name="rq3_detection_skill", dir=FUSION_DIR, file="detection_skill_bootstrap.csv",
        keys=[:product, :method, :threshold, :metric], point=:delta, lo=:ci_low, hi=:ci_high,
        p=:pvalue, verdict=row -> (excludes_zero(row.lo, row.hi),), cited=row -> false),
    (; name="fy4b_quality", dir=QUALITY_DIR, file="paired.csv",
        keys=[:scheme, :baseline, :method], point=:rel_RMSE, lo=:rel_RMSE_lo, hi=:rel_RMSE_hi,
        p=nothing, verdict=row -> (excludes_zero(row.lo, row.hi),),
        cited=row -> (row.scheme, row.baseline, row.method) in CITED_QUALITY),
]

"""Error unless the `_block1d` copy of every family file equals the published one byte for byte."""
function check_guard()
    files = vcat([(f.dir, f.file) for f in FAMILIES], [(CLAIMS_DIR, "paired_comparisons_by_method.csv")])
    for (dir, file) in unique(files)
        published, guard = joinpath(dir, file), joinpath(block_dir(dir, 1), file)
        isfile(guard) || error("no $guard: run the scripts with --block-days 1 first")
        read(published) == read(guard) || error("$guard differs from $published")
    end
    println("Guard: every --block-days 1 file equals the published one")
end

"""A family's rows at one block length, as `key => (; point, lo, hi, p)`."""
function family_values(family, dir)
    table = CSV.read(joinpath(dir, family.file), DataFrame; stringtype=String)
    result = Dict{Tuple,NamedTuple}()
    for row in eachrow(table)
        key = Tuple(string(row[k]) for k in family.keys)
        haskey(result, key) && error("$(family.file): duplicate key $key")
        result[key] = (; point=Float64(row[family.point]), lo=Float64(coalesce(row[family.lo], NaN)),
            hi=Float64(coalesce(row[family.hi], NaN)),
            p=family.p === nothing ? NaN : Float64(coalesce(row[family.p], NaN)))
    end
    return result
end

function compare_main()
    mkpath(OUTDIR)
    inputs = String[]
    for family in FAMILIES, dir in vcat(family.dir, [block_dir(family.dir, n) for n in vcat(1, BLOCKS)])
        push!(inputs, joinpath(dir, family.file))
    end
    provenance = start_provenance(ROOT, OUTDIR; inputs=unique(inputs))
    check_guard()

    rows = NamedTuple[]
    for family in FAMILIES
        day = family_values(family, family.dir)
        blocks = Dict(n => family_values(family, block_dir(family.dir, n)) for n in BLOCKS)
        for (key, d) in sort(collect(day); by=first)
            keyrow = NamedTuple{Tuple(family.keys)}(key)
            cited = family.cited(keyrow)
            verdict(v) = family.verdict(merge(keyrow, v))
            day_verdict = verdict(d)
            record = Pair{Symbol,Any}[:family => family.name, :key => join(key, " | "), :cited => cited,
                :point => d.point, :lo_day => d.lo, :hi_day => d.hi, :p_day => d.p,
                :verdict_day => join(day_verdict, "/")]
            for n in BLOCKS
                b = get(blocks[n], key, nothing)
                b === nothing && error("$(family.file) at $n days has no row $key")
                # `compare_days` sums per-block totals, so its point estimates differ from the
                # day-block run's in summation order alone (up to ~1e-12 relative).
                (isapprox(b.point, d.point; rtol=1e-9) || (isnan(b.point) && isnan(d.point))) ||
                    error("$(family.file) $key: point estimate moved at $n days")
                width_day, width = d.hi - d.lo, b.hi - b.lo
                append!(record, [Symbol("lo_$(n)d") => b.lo, Symbol("hi_$(n)d") => b.hi,
                    Symbol("p_$(n)d") => b.p, Symbol("verdict_$(n)d") => join(verdict(b), "/"),
                    Symbol("flip_$(n)d") => verdict(b) != day_verdict,
                    Symbol("widening_$(n)d") => width_day > 0 ? width / width_day : NaN])
            end
            push!(rows, (; record...))
        end
    end
    verdicts = DataFrame(rows)
    summary = combine(groupby(verdicts, :family),
        nrow => :rows, :cited => sum => :cited_rows,
        [Symbol("widening_$(n)d") => (w -> median(filter(isfinite, w))) => Symbol("median_widening_$(n)d")
            for n in BLOCKS]...,
        [Symbol("flip_$(n)d") => sum => Symbol("flips_$(n)d") for n in BLOCKS]...,
        [[Symbol("flip_$(n)d"), :cited] => ((f, c) -> sum(f .& c)) => Symbol("cited_flips_$(n)d")
            for n in BLOCKS]...)
    write_csv_atomic(joinpath(OUTDIR, "verdicts.csv"), verdicts)
    write_csv_atomic(joinpath(OUTDIR, "summary.csv"), summary)
    CSV.write(joinpath(OUTDIR, "provenance.csv"), finish_provenance!(provenance, ROOT))

    println("\nPer family (median widening of the interval; verdict flips, all and cited):")
    for row in eachrow(summary)
        @printf("  %-26s rows %3d (cited %2d)  widening %.2fx at 7 d, %.2fx at 14 d  flips %3d / %3d (cited %d / %d)\n",
            row.family, row.rows, row.cited_rows, row.median_widening_7d, row.median_widening_14d,
            row.flips_7d, row.flips_14d, row.cited_flips_7d, row.cited_flips_14d)
    end
    flipped = filter(row -> row[Symbol("flip_$(JUDGED)d")] || row.flip_14d, verdicts)
    println("\nRows whose verdict changes (cited first):")
    for row in eachrow(sort(flipped, [order(:cited; rev=true), :family]))
        @printf("  %s %-26s %-55s point %+.4f  day %+.4f..%+.4f p %.3f [%s]  7 d %+.4f..%+.4f p %.3f [%s]  14 d %+.4f..%+.4f [%s]\n",
            row.cited ? "CITED" : "     ", row.family, row.key, row.point, row.lo_day, row.hi_day, row.p_day,
            row.verdict_day, row.lo_7d, row.hi_7d, row.p_7d, row.verdict_7d, row.lo_14d, row.hi_14d,
            row.verdict_14d)
    end
    cited_flips = count(row -> row.cited && row[Symbol("flip_$(JUDGED)d")], eachrow(verdicts))
    println(cited_flips == 0 ?
        "\nVerdict: no cited verdict flips at $JUDGED days; the paper's verdicts stand" :
        "\nVerdict: $cited_flips cited verdict(s) flip at $JUDGED days; the user decides the wording")
    println("Wrote $(relpath(OUTDIR, ROOT))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_main()
end
