#!/usr/bin/env julia

#=
Do calibrated (`GPM`, IMERG Final) and uncalibrated (`GPM_uncal`) IMERG read the same pixel at
every gauge?

They come from different extractions. `GPM` was sampled in Earth Engine (`reduceRegions`, scale
11132, on the asset's Float32 coordinates); `GPM_uncal` is rebuilt from the GES DISC grid by
`IMERGUncalStations` (`floor((v - c0) / 0.1 + 0.5)` on `station_meta.csv`). For a gauge on a 0.1 deg
cell edge the two rules can fall on different sides.

Earth Engine's pixel cannot be read directly - there is no calibrated grid on disk - but it can be
identified. IMERG Final is IMERG Uncal rescaled by a monthly gauge adjustment, so at the matching
pixel the two series share their dry hours and keep a near-constant ratio within a month. For each
study gauge this scores the 3x3 neighbourhood of Uncal pixels against the calibrated series and
takes the best one; the interior gauges, where every rule agrees, validate the method.

It then recomputes the storm-core calibration ratio (`final_over_uncal`, cited as rq1_cores `cal`)
with the Uncal value read at Earth Engine's pixel. The assertions at the end hold the numbers the
caveat on `prep_gpm_uncal_stations` in paper/sources.toml quotes; if one fails, update that note.

    julia --project=. scripts/verify_imerg_pixel_agreement.jl

Reads the raw grid (about 25 s) and writes CSVs to output/imerg_pixel_agreement/.
=#

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("IMERGUncalStations", "TableIO", "RunProvenance")

using CSV, DataFrames, Dates, NCDatasets, Printf, Statistics
using Main.IMERGUncalStations: imerg_cell, imerg_halfhour_files
using Main.TableIO: write_csv_atomic
using Main.RunProvenance: start_provenance, finish_provenance!

const RAW_DIR = joinpath(ROOT, "data", "raw", "gpm_imerg_v07_uncal_hubei_2022_2024")
const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const CALIBRATED = joinpath(STUDY_DATA, "hubei_gpm_hourly_2022_2024_full_aligned.csv")
const STATION_META = joinpath(STUDY_DATA, "station_meta.csv")
const CORE_GAUGES = joinpath(ROOT, "output", "heavy_core_diagnostics_gpmuncal", "core_gauges.csv")
const OUTDIR = joinpath(ROOT, "output", "imerg_pixel_agreement")
const DAYS = Date(2022, 1, 1):Day(1):Date(2024, 12, 31)

"""Hourly Uncal over the whole subset grid, `cube[lat, lon, hour]`; hour `k` starts at UTC 2022-01-01 + k-1 h."""
function hourly_cube()
    lon, lat = NCDataset(last(first(imerg_halfhour_files(RAW_DIR, first(DAYS))))) do ds
        Float64.(ds["lon"][:]), Float64.(ds["lat"][:])
    end
    read_grid(path) = NCDataset(ds -> Float64.(coalesce.(ds["precipitationUncal"][:, :, 1], NaN32)), path)
    cube = Array{Float32}(undef, length(lat), length(lon), 24 * length(DAYS))
    k = 0
    for day in DAYS
        files = imerg_halfhour_files(RAW_DIR, day)
        for h in 1:2:47
            k += 1
            cube[:, :, k] = Float32.((read_grid(files[h][2]) .+ read_grid(files[h+1][2])) ./ 2)
        end
    end
    return lon, lat, cube
end

on_edge(v) = round(Int, v * 1000) % 100 == 0

"""Dry/wet agreement, wet-hour correlation and median within-month spread of log(calibrated / uncalibrated)."""
function agreement(c, u, months)
    ok = isfinite.(c) .& isfinite.(u)
    c, u, months = c[ok], u[ok], months[ok]
    wet = (c .> 0) .| (u .> 0)
    spread = Float64[]
    for month in unique(months)
        both = (months .== Ref(month)) .& (c .> 0) .& (u .> 0)
        count(both) >= 10 || continue
        ratio = log.(c[both] ./ u[both])
        push!(spread, median(abs.(ratio .- median(ratio))))
    end
    return (; dry_wet_agreement=mean((c .> 0) .== (u .> 0)), wet_hour_r=cor(c[wet], u[wet]),
        ratio_spread=isempty(spread) ? NaN : median(spread))
end

function main()
    provenance = start_provenance(ROOT, OUTDIR; inputs=[CALIBRATED, STATION_META, CORE_GAUGES])
    lon, lat, cube = hourly_cube()
    lon0, lat0 = imerg_cell(lon[1], -179.95), imerg_cell(lat[1], -89.95)
    cell(lonv, latv) = (imerg_cell(latv, -89.95) - lat0 + 1, imerg_cell(lonv, -179.95) - lon0 + 1)

    calibrated = CSV.read(CALIBRATED, DataFrame; types=Dict(:time => String))
    ids = names(calibrated)[2:end]
    meta = CSV.read(STATION_META, DataFrame; types=Dict(:station_id => String))
    meta = meta[indexin(ids, meta.station_id), :]
    # Rows are labelled UTC hour start + 9 h, like the Uncal file.
    utc = [DateTime(chop(t)) - Hour(9) for t in calibrated.time]
    hour_index = [Int(Dates.value(t - DateTime(first(DAYS))) ÷ 3_600_000) + 1 for t in utc]
    months = yearmonth.(utc)

    scores = NamedTuple[]
    for (k, id) in enumerate(ids)
        c = Float64.(coalesce.(calibrated[!, id], NaN))
        i0, j0 = cell(meta.lon[k], meta.lat[k])
        for di in -1:1, dj in -1:1
            i, j = i0 + di, j0 + dj
            (1 <= i <= length(lat) && 1 <= j <= length(lon)) || continue
            push!(scores, (; station_id=id, lon=meta.lon[k], lat=meta.lat[k], edge_lon=on_edge(meta.lon[k]),
                edge_lat=on_edge(meta.lat[k]), di, dj, agreement(c, Float64.(cube[i, j, hour_index]), months)...))
        end
    end
    scores = DataFrame(scores)
    write_csv_atomic(joinpath(OUTDIR, "pixel_scores.csv"), scores)

    # Earth Engine's pixel: the best dry/wet agreement, then correlation. (di, dj) is its offset
    # from the pixel GPM_uncal reads.
    gauges = combine(groupby(scores, :station_id)) do g
        order = sortperm(collect(zip(g.dry_wet_agreement, g.wet_hour_r)); rev=true)
        best, own = g[order[1], :], only(eachrow(g[(g.di .== 0) .& (g.dj .== 0), :]))
        (; lon=best.lon, lat=best.lat, edge=best.edge_lon || best.edge_lat, di=best.di, dj=best.dj,
           spread_earth_engine_pixel=best.ratio_spread, spread_uncal_pixel=own.ratio_spread,
           agreement_earth_engine_pixel=best.dry_wet_agreement, agreement_uncal_pixel=own.dry_wet_agreement)
    end
    gauges.same_pixel = (gauges.di .== 0) .& (gauges.dj .== 0)
    write_csv_atomic(joinpath(OUTDIR, "gauge_pixels.csv"), gauges)

    interior, edge = gauges[.!gauges.edge, :], gauges[gauges.edge, :]
    mismatched = edge[.!edge.same_pixel, :]
    @printf("Interior gauges: %d, own pixel best for %d (ratio spread median %.3f, max %.3f)\n", nrow(interior),
        count(interior.same_pixel), median(interior.spread_earth_engine_pixel), maximum(interior.spread_earth_engine_pixel))
    @printf("Edge gauges: %d; same pixel %d, neighbouring pixel %d\n", nrow(edge), count(edge.same_pixel), nrow(mismatched))
    @printf("  at the %d: ratio spread median %.3f at Earth Engine's pixel (max %.3f) vs %.3f at GPM_uncal's (min %.3f)\n",
        nrow(mismatched), median(mismatched.spread_earth_engine_pixel), maximum(mismatched.spread_earth_engine_pixel),
        median(mismatched.spread_uncal_pixel), minimum(mismatched.spread_uncal_pixel))
    offsets = sort(combine(groupby(mismatched, [:di, :dj]), nrow => :gauges), :gauges; rev=true)
    println(offsets)

    # The storm-core calibration ratio with Uncal at Earth Engine's pixel. The 08-08 BJT event day
    # is UTC day `event_day`.
    offset = Dict(r.station_id => (r.di, r.dj) for r in eachrow(gauges))
    core = CSV.read(CORE_GAUGES, DataFrame; types=Dict(:station_id => String))
    function day_total(r, di, dj)
        i, j = cell(r.lon, r.lat)
        d = Dates.value(Date(r.event_day) - first(DAYS)) + 1
        return sum(Float64.(cube[i+di, j+dj, 24(d-1)+1:24d]))
    end
    own = [day_total(r, 0, 0) for r in eachrow(core)]
    maximum(abs.(own .- core.GPM_uncal_mm)) < 1e-4 ||
        error("the rebuilt event totals do not reproduce core_gauges.csv's GPM_uncal_mm")
    core.uncal_at_earth_engine_pixel = [day_total(r, offset[r.station_id]...) for r in eachrow(core)]
    core.mismatched = [offset[id] != (0, 0) for id in core.station_id]
    ratio = combine(groupby(core, [:event_day, :event_type]),
        nrow => :n_core, :mismatched => count => :n_mismatched,
        [:GPM_mm, :GPM_uncal_mm] => ((g, u) -> mean(g) / mean(u)) => :final_over_uncal,
        [:GPM_mm, :uncal_at_earth_engine_pixel] => ((g, u) -> mean(g) / mean(u)) => :final_over_uncal_earth_engine_pixel)
    write_csv_atomic(joinpath(OUTDIR, "calibration_ratio.csv"), ratio)
    println(ratio)

    # The method is only trusted if the interior gauges validate it and the two pixels separate.
    all(interior.same_pixel) || error("an interior gauge's own pixel does not score best; the method cannot tell")
    maximum(mismatched.spread_earth_engine_pixel) < minimum(mismatched.spread_uncal_pixel) ||
        error("the two candidate pixels' ratio spreads overlap; the mismatches are not decisive")
    # The numbers the caveat in paper/sources.toml quotes. If one moves, update that note.
    (nrow(gauges), nrow(edge), nrow(mismatched)) == (237, 73, 49) || error("gauge counts moved")
    Dict(Tuple(r[[:di, :dj]]) => r.gauges for r in eachrow(offsets)) == Dict((-1, 0) => 33, (0, -1) => 13, (-1, -1) => 3) ||
        error("the mismatch directions moved")
    row = only(eachrow(ratio[string.(ratio.event_day) .== "2024-07-30", :]))
    (round(row.final_over_uncal; digits=3), round(row.final_over_uncal_earth_engine_pixel; digits=3)) == (1.681, 1.724) ||
        error("the cited 2024-07-30 calibration ratio moved")

    write_csv_atomic(joinpath(OUTDIR, "provenance.csv"), finish_provenance!(provenance, ROOT))
    println("Every check holds; tables in $OUTDIR")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
