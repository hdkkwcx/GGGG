#!/usr/bin/env julia
module FY4BPreprocessing

using NCDatasets, CSV, DataFrames, Dates
using Main.TableIO: write_csv_atomic

export aggregate_fy4b_hourly, aggregate_fy4b_variants, build_hourly_qc, find_nc_files,
    parse_nc_filename

# =========================
# FY4B Aggregation: 15min NetCDF -> Hourly CSV
# =========================

# Two levels up, because this file lives in `src/sources/`. The consts below are defaults the
# `prepare_fy4b_*` scripts override; they are string consts, so a wrong root here would not fail
# on load - it would surface much later, as a missing input inside whichever script read through
# one of them.
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const FY4B_DATA_DIR = joinpath(ROOT, "data", "FY4B")
const STATION_META_FILE = joinpath(ROOT, "data", "hubei_station_meta.csv")
const OUTPUT_FILE = joinpath(ROOT, "data", "processed", "hubei_fy4b_hourly_2022_2025_JunSep_strict_navcorrected.csv")
const QC_FILE = joinpath(ROOT, "output", "input_audit", "fy4b_hourly_qc_2022_2025_JunSep_strict_navcorrected.csv")

# Satellite parameters for FY4B (from NC files)
const DEFAULT_SAT_LON = 133.0   # Fallback satellite longitude (degrees)
const SAT_HEIGHT = 42164.0e3    # Satellite height from Earth center (m)
const EARTH_RADIUS = 6378.137e3 # Earth equatorial radius (m)
const GRID_SIZE = 2748          # Grid dimensions
const RESOLUTION = 4000.0       # 4km at nadir (m)
const STATION_TIME_OFFSET = Hour(8) # FY4B filenames are UTC; station timestamps are Beijing time

# Computed constants
const SAT_DISTANCE = SAT_HEIGHT - EARTH_RADIUS  # Distance from surface to satellite

# The module's code, split by concern in `fy4b/`:
# - `fy4b_projection.jl`, `fy4b_filenames.jl`: from a lon/lat to a disk pixel, and what a file
#   name says (its 15 minutes, its satellite longitude).
# - `fy4b_station_extract.jl`, `fy4b_hourly.jl`: the single-pixel series and the pipeline that
#   writes the shipped hourly file (`aggregate_fy4b_hourly`).
# - `fy4b_variants.jl`, `fy4b_offset_grid.jl`, `fy4b_rebuild.jl`: the rebuild in window variants,
#   its grid of pixel shifts, and the pipeline that writes them (`aggregate_fy4b_variants`).
# - `fy4b_files_qc.jl`: which raw files there are, which open, and which hours are complete.
include("fy4b/fy4b_projection.jl")
include("fy4b/fy4b_filenames.jl")
include("fy4b/fy4b_station_extract.jl")
include("fy4b/fy4b_variants.jl")
include("fy4b/fy4b_offset_grid.jl")
include("fy4b/fy4b_rebuild.jl")
include("fy4b/fy4b_files_qc.jl")
include("fy4b/fy4b_hourly.jl")

end # module
