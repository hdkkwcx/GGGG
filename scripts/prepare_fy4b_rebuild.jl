#!/usr/bin/env julia

# Rebuild the FY4B station series from the raw 15-minute QPE files, in several gauge-free
# variants from one pass over the files (`FY4BPreprocessing.aggregate_fy4b_variants`).
#
#   julia --project=. scripts/prepare_fy4b_rebuild.jl                  # 2022-2024 -> output/fy4b_rebuild_2022_2024/
#   julia --project=. scripts/prepare_fy4b_rebuild.jl --month 2022-06  # one month -> output/fy4b_rebuild_test_202206/
#
# Variants are named w<window>_h<cloud-top km>. `w1_h0` is the shipped method - nearest pixel, no
# parallax - and must reproduce
# data/processed/study_area/hubei_fy4b_hourly_2022_2024_full_strict_navcorrected.csv on its hours.
# `w3`/`w5` average a 3x3/5x5 pixel box (about 12/20 km). `w3_h10` samples where a 10 km cloud top
# over the gauge appears in the image; it is kept only as a reported check, because on June 2022
# it scored worse than no shift (r 0.343 vs 0.351) and the displacement the gauges show is not
# parallax-shaped (2022 west, 2024 east). All are fixed in advance, none tuned on gauges, so the
# series can enter the benchmark's cross-validation without leaking.
#
# fy4b_offset_grid.nc holds the 3x3 mean at every pixel shift up to +-6: the candidates from which
# the benchmark picks a shift per fold on training gauges only. Hours stay strict: all four
# 15-minute files, or no value. The 1-, 3- and 24-hour files in the raw tree are gauge-corrected
# and never read.
#
# The raw tree (data/FY4B, a link to the external drive) is only read. Like the IMERG grid it is
# not hashed; every file used is listed in fy4b_hourly_qc.csv.

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("FY4BPreprocessing", "RunProvenance")

using CSV, Dates
using Main.FY4BPreprocessing
using Main.RunProvenance: start_provenance, finish_provenance!

const RAW_DIR = joinpath(ROOT, "data", "FY4B")
const STATION_META = joinpath(ROOT, "data", "processed", "study_area", "station_meta.csv")
const VARIANTS = [
    (name="w1_h0", window=1, height=0.0),
    (name="w3_h0", window=3, height=0.0),
    (name="w5_h0", window=5, height=0.0),
    (name="w3_h10", window=3, height=10_000.0),
]
const OFFSET_GRID = (window=3, radius=6)

function main(args=ARGS)
    index = findfirst(==("--month"), args)
    if index === nothing
        years, months = 2022:2024, 1:12
        outdir = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024")
    else
        month_start = Date(args[index + 1], dateformat"yyyy-mm")
        years, months = year(month_start):year(month_start), [month(month_start)]
        outdir = joinpath(ROOT, "output", "fy4b_rebuild_test_" * Dates.format(month_start, "yyyymm"))
    end
    isdir(RAW_DIR) || error("Raw FY4B tree not found at $RAW_DIR - is the external drive mounted?")
    provenance = start_provenance(ROOT, outdir; inputs=[STATION_META])
    started = time()
    aggregate_fy4b_variants(;
        data_dir=RAW_DIR, station_meta_file=STATION_META, outdir, years, months,
        variants=VARIANTS, offset_grid=OFFSET_GRID,
    )
    println("Rebuild took $(round((time() - started) / 60; digits=1)) min")
    CSV.write(joinpath(outdir, "provenance.csv"), finish_provenance!(provenance, ROOT))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
