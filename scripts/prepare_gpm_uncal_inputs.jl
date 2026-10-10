#!/usr/bin/env julia

# Build the uncalibrated IMERG (V07 `precipitationUncal`) benchmark input the same way the
# calibrated `hubei_gpm_hourly_2022_2024_full_aligned.csv` is built: the monthly station series
# are shifted +9 h and widened by `prepare_satellite_wide`, then cut to the study-area stations.
# Read by `run_interpolation_benchmark.jl --gpm-uncal`.

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("MGERDataPrep", "StudyArea")

using CSV, DataFrames
using Main.MGERDataPrep
using Main.StudyArea

const PROCESSED = joinpath(ROOT, "data", "processed")
const STUDY_DATA = joinpath(PROCESSED, "study_area")
const UNCAL_DIR = joinpath(ROOT, "output", "gpm_imerg_uncal_stations_2022_2024")
const FILE_NAME = "hubei_gpm_uncal_hourly_2022_2024_full_aligned.csv"

function main()
    mkpath(joinpath(ROOT, "output", "input_audit"))
    wide = prepare_satellite_wide(
        "gpm", UNCAL_DIR, joinpath(PROCESSED, FILE_NAME);
        value_column=:gpm_mm_h,
        years=2022:2024,
        months=1:12,
        duplicate_qc_file=joinpath(ROOT, "output", "input_audit", "gpm_uncal_full_duplicate_station_times.csv"),
    )
    println(wide)
    stations = filter_stations(CSV.read(joinpath(ROOT, "data", "hubei_station_meta.csv"), DataFrame))
    study = StudyArea._subset_wide_file(
        joinpath(PROCESSED, FILE_NAME), joinpath(STUDY_DATA, FILE_NAME), Vector{String}(stations.station_id),
    )
    println(study)
    return (; wide, study)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
