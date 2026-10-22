#!/usr/bin/env julia

# Build the monthly uncalibrated IMERG (V07 `precipitationUncal`) station series that
# `prepare_gpm_uncal_inputs.jl` and `run_no_rain_evaluation.jl` read, from the half-hourly GES DISC
# subsets in data/raw/gpm_imerg_v07_uncal_hubei_2022_2024/ (`IMERGUncalStations`).
#
#   julia --project=. scripts/prepare_gpm_uncal_stations.jl             # output/gpm_imerg_uncal_stations_2022_2024/
#   julia --project=. scripts/prepare_gpm_uncal_stations.jl <outdir>
#
# One `gpm_hubei_hourly_long_<yyyymm>.csv` per month of 2022-2024, every gauge of
# data/hubei_station_meta.csv that lies on the subset grid, plus provenance.csv. Widened by
# `prepare_gpm_uncal_inputs.jl`, the series reproduces
# data/processed/study_area/hubei_gpm_uncal_hourly_2022_2024_full_aligned.csv byte for byte
# (checked 2026-10-22).

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("IMERGUncalStations", "RunProvenance")

using CSV, DataFrames
using Main.IMERGUncalStations
using Main.RunProvenance: start_provenance, finish_provenance!

const RAW_DIR = joinpath(ROOT, "data", "raw", "gpm_imerg_v07_uncal_hubei_2022_2024")
const STATION_META = joinpath(ROOT, "data", "hubei_station_meta.csv")
const DEFAULT_OUTDIR = joinpath(ROOT, "output", "gpm_imerg_uncal_stations_2022_2024")

function main(args=ARGS)
    outdir = isempty(args) ? DEFAULT_OUTDIR : abspath(args[1])
    # The raw grid (52,608 files) is not hashed: its folder and request are documented in
    # paper/sources.toml, and the series is checked against the prepared file it feeds.
    provenance = start_provenance(ROOT, outdir; inputs=[STATION_META])
    stations = CSV.read(STATION_META, DataFrame; types=Dict(:station_id => String))
    outside = write_imerg_uncal_months(RAW_DIR, stations, outdir)
    isempty(outside) || println("Left out, off the subset grid: $(join(outside, ", "))")
    CSV.write(joinpath(outdir, "provenance.csv"), finish_provenance!(provenance, ROOT))
    println("Wrote $(nrow(stations) - length(outside)) stations x 36 months to $outdir")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
