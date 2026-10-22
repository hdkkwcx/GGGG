#!/usr/bin/env julia

# Restore data/dem_ShiYan_1km.tif, the relief behind paper Fig. 1 (scripts/plot_study_area_map.py).
#
#   julia --project=. scripts/download_dem_shiyan.jl
#   julia --project=. scripts/download_dem_shiyan.jl <output.tif>
#
# scripts/clip_dem.R made this file from a nationwide 1 km SRTM DEM that is not in this project
# (it read it from the VICTools R package's database on a local drive). The clipped file itself
# is published in the upstream repository jl-pkgs/GeoWeightedRegression.jl; it is fetched from a
# pinned commit and accepted only if its SHA-256 matches the one Fig. 1's plot_provenance.csv
# records. An existing file is never overwritten: it is checked against the same hash instead.

using Downloads
using SHA

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_OUTPUT = joinpath(ROOT, "data", "dem_ShiYan_1km.tif")
const URL = "https://raw.githubusercontent.com/jl-pkgs/GeoWeightedRegression.jl/" *
    "d277481c03d0c1dca18309d966d6cd9e8d1a94f3/data/dem_ShiYan_1km.tif"
const EXPECTED_SHA256 = "7e7d6587a176bab3ed626985228da73efb3067b2d20b8c53c54062c3a223cebc"

file_sha256(path) = open(io -> bytes2hex(sha256(io)), path)

function main(args=ARGS)
    output = isempty(args) ? DEFAULT_OUTPUT : abspath(args[1])
    if isfile(output)
        file_sha256(output) == EXPECTED_SHA256 ||
            error("$output exists but its SHA-256 differs from the one Fig. 1 was drawn from; " *
                  "move it aside to restore the published copy")
        println("$output is already the published DEM (SHA-256 matches)")
        return
    end
    tmp = Downloads.download(URL)
    got = file_sha256(tmp)
    if got != EXPECTED_SHA256
        rm(tmp)
        error("downloaded file has SHA-256 $got, expected $EXPECTED_SHA256")
    end
    mkpath(dirname(output))
    mv(tmp, output)
    println("wrote $output (SHA-256 verified)")
end

main()
