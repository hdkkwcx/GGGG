"""
Hourly IMERG V07 `precipitationUncal` at the gauges, from the half-hourly GES DISC subsets.

The paper's uncalibrated IMERG input (`GPM_uncal`) starts as monthly long tables of this series,
which `prepare_gpm_uncal_inputs.jl` widens exactly as the calibrated Earth Engine exports are
widened. This module rebuilds those tables from `data/raw/gpm_imerg_v07_uncal_hubei_2022_2024/`,
the 0.1 deg, 30-min subsets (`<year>/<day of year>/3B-HHR.MS.MRG.3IMERG.<yyyymmdd>-S<hhmmss>-...nc4`).

Two choices are not free: they are the ones the original series was built with, recovered by
matching it value for value at all 237 study gauges.

- A gauge reads the cell `floor((v - c0) / 0.1 + 0.5)` along each axis, where `c0` is the global
  grid's first cell centre (-179.95 E, -89.95 N), evaluated in Float64 on the station table's
  coordinates. It is the nearest cell centre; for a gauge on a cell edge, the floating-point
  rounding of this exact expression decides the side, so it must not be rewritten.
- An hour (UTC, labelled by its start) is the mean of its two half-hourly rates, each widened to
  Float64 before averaging. A missing half (fill value) makes the hour NaN.
"""
module IMERGUncalStations

using NCDatasets, DataFrames, Dates
using Main.TableIO: write_csv_atomic

export imerg_cell, imerg_halfhour_files, imerg_hourly_station_series, write_imerg_uncal_months

const RESOLUTION = 0.1
const FIRST_CENTRE_LON = -179.95
const FIRST_CENTRE_LAT = -89.95

"""Global 0-based index of the cell `v` reads along an axis whose first centre is `first`."""
imerg_cell(v::Real, first::Real) = floor(Int, (Float64(v) - first) / RESOLUTION + 0.5)

"""
`(time, path)` of the 48 half-hourly files of UTC day `day` under `dir`, sorted by time.

Times come from the file names: the files' own `time` variable declares a Julian calendar and
decodes 13 days late.
"""
function imerg_halfhour_files(dir::AbstractString, day::Date)
    daydir = joinpath(dir, string(year(day)), lpad(dayofyear(day), 3, '0'))
    stamp = Dates.format(day, "yyyymmdd")
    files = Tuple{DateTime,String}[]
    for file in readdir(daydir)
        m = match(r"\.(\d{8})-S(\d{6})-.*\.nc4$", file)
        (m === nothing || m.captures[1] != stamp) && continue
        push!(files, (DateTime(m.captures[1] * m.captures[2], dateformat"yyyymmddHHMMSS"),
            joinpath(daydir, file)))
    end
    length(files) == 48 || error("expected 48 half-hourly IMERG files for $day in $daydir, found $(length(files))")
    return sort!(files)
end

"""The `precipitationUncal` rate (mm/h) at `cells` (`(lat, lon)` indices), NaN where filled."""
function read_rates(path::AbstractString, cells)
    grid = NCDataset(path) do ds
        dimnames(ds["precipitationUncal"]) == ("lat", "lon", "time") || error("unexpected dimensions in $path")
        ds["precipitationUncal"][:, :, 1]
    end
    return [Float64(coalesce(grid[i, j], NaN32)) for (i, j) in cells]
end

"""
The month `(year, month)` of hourly station series: a long table of `time` (UTC hour start,
`yyyy-mm-ddTHH:MM:SSZ`), `station_id` and `gpm_mm_h`, in the shape `prepare_satellite_wide`
reads. `stations` needs `station_id`, `lon` and `lat`. Returns `(; table, outside)`, `outside`
being the stations that fall off the subset grid and are left out.
"""
function imerg_hourly_station_series(dir::AbstractString, stations::AbstractDataFrame, year::Int, month::Int)
    days = Date(year, month, 1):Day(1):lastdayofmonth(Date(year, month, 1))
    lon, lat = NCDataset(last(first(imerg_halfhour_files(dir, first(days))))) do ds
        Float64.(ds["lon"][:]), Float64.(ds["lat"][:])
    end
    lon0, lat0 = imerg_cell(lon[1], FIRST_CENTRE_LON), imerg_cell(lat[1], FIRST_CENTRE_LAT)
    ids, cells, outside = String[], Tuple{Int,Int}[], String[]
    for row in eachrow(stations)
        i = imerg_cell(row.lat, FIRST_CENTRE_LAT) - lat0 + 1
        j = imerg_cell(row.lon, FIRST_CENTRE_LON) - lon0 + 1
        if 1 <= i <= length(lat) && 1 <= j <= length(lon)
            push!(ids, string(row.station_id)); push!(cells, (i, j))
        else
            push!(outside, string(row.station_id))
        end
    end
    times, station_ids, values = String[], String[], Float64[]
    for day in days
        files = imerg_halfhour_files(dir, day)
        for k in 1:2:47
            (start, first_half), (next, second_half) = files[k], files[k+1]
            minute(start) == 0 && next == start + Minute(30) || error("unpaired half hours at $start")
            hour_value = (read_rates(first_half, cells) .+ read_rates(second_half, cells)) ./ 2
            append!(times, fill(Dates.format(start, "yyyy-mm-ddTHH:MM:SS") * "Z", length(ids)))
            append!(station_ids, ids)
            append!(values, hour_value)
        end
    end
    return (; table=DataFrame(time=times, station_id=station_ids, gpm_mm_h=values), outside)
end

"""
Write one `gpm_hubei_hourly_long_<yyyymm>.csv` per month of `years` into `outdir`, the file names
`prepare_gpm_uncal_inputs.jl` reads. Returns the stations left out as off the grid.
"""
function write_imerg_uncal_months(dir::AbstractString, stations::AbstractDataFrame, outdir::AbstractString;
                                  years=2022:2024)
    mkpath(outdir)
    outside = Set{String}()
    for y in years, m in 1:12
        series = imerg_hourly_station_series(dir, stations, y, m)
        union!(outside, series.outside)
        write_csv_atomic(joinpath(outdir, "gpm_hubei_hourly_long_$(y)$(lpad(m, 2, '0')).csv"), series.table)
        println("  $(y)-$(lpad(m, 2, '0')): $(nrow(series.table)) rows")
    end
    return sort!(collect(outside))
end

end # module
