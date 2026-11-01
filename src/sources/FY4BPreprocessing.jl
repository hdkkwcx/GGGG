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

"""
    latlon_to_scan_angles(lat, lon, sat_lon=DEFAULT_SAT_LON; height=0.0)

Convert latitude/longitude to scan angles (x, y in radians) for geostationary projection.
Uses standard geostationary projection formula (similar to GOES-R ABI).
Returns: (scan_x, scan_y) in radians, or (NaN, NaN) if not visible.

`height` (metres) lifts the point along the local vertical before projecting, which gives where
a cloud top at that height over the spot appears in the image: the parallax-displaced position.
Over Hubei a 10 km top is displaced about 10 km toward the NW from 133E and about 8 km toward
the N from 105E. `height == 0` is the surface point, computed exactly as before.
"""
function latlon_to_scan_angles(lat, lon, sat_lon=DEFAULT_SAT_LON; height=0.0)
    # Convert to radians
    φ = deg2rad(lat)
    λ = deg2rad(lon - sat_lon)  # Relative longitude

    # WGS84 parameters
    re = EARTH_RADIUS           # Equatorial radius
    rp = 6356.75231414e3        # Polar radius (m)
    f = (re - rp) / re          # Flattening
    e2 = f * (2 - f)            # Eccentricity squared

    # Geocentric latitude (accounting for Earth flattening)
    φ_c = atan((1 - e2) * tan(φ))

    # Distance from Earth center to surface point
    rc = rp / sqrt(1 - e2 * cos(φ_c)^2)

    # Cartesian coordinates of surface point (Earth-centered, x-axis through 0°E)
    xs = rc * cos(φ_c) * cos(λ)
    ys = rc * cos(φ_c) * sin(λ)
    zs = rc * sin(φ_c)

    # Lift along the geodetic normal; skipped at zero so the surface path stays bit-identical.
    if height != 0
        xs += height * cos(φ) * cos(λ)
        ys += height * cos(φ) * sin(λ)
        zs += height * sin(φ)
    end

    # Visibility check: point must be on Earth side facing the satellite
    # The satellite is at (SAT_HEIGHT, 0, 0) in this coordinate system
    # Point is visible if it's not behind the Earth's limb
    dot_product = xs * (xs - SAT_HEIGHT) + ys^2 + zs^2
    if dot_product > 0
        # Check if behind Earth
        # Actually, we need to check if the line from satellite to point intersects Earth
        # Simplified: check the angle from satellite
        d2 = (SAT_HEIGHT - xs)^2 + ys^2 + zs^2
        cos_angle = (SAT_HEIGHT - xs) / sqrt(d2)
        if cos_angle < 0.156  # cos(81°) - maximum view angle
            return (NaN, NaN)
        end
    end

    # Scan angles from satellite perspective
    # x (E-W, west negative): angle in the equatorial plane
    scan_x = atan(ys / (SAT_HEIGHT - xs))

    # y (N-S, north positive): angle perpendicular to equatorial plane
    scan_y = atan(zs / sqrt((SAT_HEIGHT - xs)^2 + ys^2))

    return (scan_x, scan_y)
end

"""
    scan_angles_to_grid(scan_x, scan_y)

Convert scan angles (radians) to FY4B grid indices (0-based).
The grid is centered at the satellite sub-point (nadir).
"""
function scan_angles_to_grid(scan_x, scan_y)
    if isnan(scan_x) || isnan(scan_y)
        return (NaN, NaN)
    end

    # Convert scan angle to linear distance at nadir
    # At nadir: angle = distance / (H - R)
    # scan_angle = distance / SAT_DISTANCE
    # distance = scan_angle * SAT_DISTANCE

    # Grid spacing in radians (at nadir, RESOLUTION meters = scan_angle * SAT_DISTANCE)
    grid_spacing_rad = RESOLUTION / SAT_DISTANCE

    # Grid center
    center = (GRID_SIZE - 1) / 2.0

    # Convert scan angles to grid indices
    # Negative scan_x (west of nadir) -> smaller x index
    # Positive scan_y (north of nadir) -> smaller y index (image coordinates)
    x_idx = center + scan_x / grid_spacing_rad
    y_idx = center - scan_y / grid_spacing_rad

    return (x_idx, y_idx)
end

"""
    parse_satellite_lon(filename)

Parse FY4B satellite subpoint longitude from the filename.
Pattern segment: N_DISK_1330E or N_DISK_1050E.
Returns degrees east as Float64.
"""
function parse_satellite_lon(filename)
    m = match(r"N_DISK_(\d{4})E", filename)
    if m === nothing
        @warn "Could not parse satellite longitude from filename; using fallback" filename fallback=DEFAULT_SAT_LON
        return DEFAULT_SAT_LON
    end
    return parse(Int, m.captures[1]) / 10.0
end

"""
    parse_nc_filename(filename)

Parse FY4B NetCDF filename to extract start/end datetime and satellite longitude.
Pattern: FY4B-_AGRI--_N_DISK_1330E_L2-_QPE-_MULT_NOM_YYYYMMDDhhmmss_YYYYMMDDhhmmss_4000M_V0001.NC
Returns: (start_datetime, end_datetime, satellite_longitude)
"""
function parse_nc_filename(filename)
    m = match(r"NOM_(\d{14})_(\d{14})_", filename)
    if m === nothing
        return nothing
    end
    start_str, end_str = m.captures
    start_dt = DateTime(start_str, dateformat"yyyymmddHHMMSS")
    end_dt = DateTime(end_str, dateformat"yyyymmddHHMMSS")
    return (start_dt, end_dt, parse_satellite_lon(filename))
end

"""
    is_15min_file(filename)

Check if file represents a 15-minute aggregation window.
True 15min files have a time span of exactly 15 minutes (0:14:59).
"""
function is_15min_file(filename)
    dt_info = parse_nc_filename(filename)
    if dt_info === nothing
        return false
    end
    start_dt, end_dt, _ = dt_info
    return second(start_dt) == 0 &&
           minute(start_dt) in (0, 15, 30, 45) &&
           end_dt == start_dt + Minute(15) - Second(1)
end

function dataset_satellite_lon(ds, filename)
    if haskey(ds, "nominal_satellite_subpoint_lon")
        value = ds["nominal_satellite_subpoint_lon"][1]
        if !ismissing(value) && isfinite(value)
            return Float64(value)
        end
    end
    return parse_satellite_lon(filename)
end

aligned_hour_end(hour_start::DateTime) = hour_start + Hour(1) + STATION_TIME_OFFSET

"""
    latlon_to_xy(lat, lon, sat_lon=DEFAULT_SAT_LON; height=0.0)

Convert latitude/longitude to FY4B fixed grid x/y coordinates.
Uses geostationary projection formulas.
Returns: (x, y) in grid coordinates (0-based index)

`height` is the cloud-top height (metres) of the parallax correction; see `latlon_to_scan_angles`.
"""
function latlon_to_xy(lat, lon, sat_lon=DEFAULT_SAT_LON; height=0.0)
    scan_x, scan_y = latlon_to_scan_angles(lat, lon, sat_lon; height)
    return scan_angles_to_grid(scan_x, scan_y)
end

"""
    station_grid_coords(stations_df, sat_lon, coord_cache)

Compute and cache station grid coordinates for a satellite subpoint longitude.
The default aggregation keeps the 133.0E projection used by the source metadata
and the published validation run. Filename-based 105.0E/133.0E re-extraction is
handled as a separate audit candidate in aggregate_fy4b_python.py.
"""
function station_grid_coords(stations_df, sat_lon, coord_cache)
    if !haskey(coord_cache, sat_lon)
        x_coords = Float64[]
        y_coords = Float64[]
        sizehint!(x_coords, nrow(stations_df))
        sizehint!(y_coords, nrow(stations_df))

        for row in eachrow(stations_df)
            x, y = latlon_to_xy(row.lat, row.lon, sat_lon)
            push!(x_coords, x)
            push!(y_coords, y)
        end

        coord_cache[sat_lon] = (x_coords, y_coords)
    end

    return coord_cache[sat_lon]
end

"""
    load_stations(meta_file)

Load station metadata.
Returns DataFrame with station_id, lon, lat.
"""
function load_stations(meta_file)
    df = CSV.read(meta_file, DataFrame)
    # Ensure correct column names
    rename!(df, lowercase.(names(df)))
    return df
end

"""
    extract_precipitation(nc_file, stations_df, coord_cache)

Extract precipitation values from NetCDF file at station locations.
Returns: Vector of precipitation values (mm) matching stations_df order.
"""
function extract_precipitation(nc_file, stations_df, coord_cache)
    ds = Dataset(nc_file)
    try
        precip = ds["Precipitation"]
        fill_val = Float32(65534.0)  # From NC attributes
        quality = haskey(ds, "DQF") ? ds["DQF"] : nothing

        sat_lon = dataset_satellite_lon(ds, basename(nc_file))
        x_coords, y_coords = station_grid_coords(stations_df, sat_lon, coord_cache)

        values = Float64[]
        sizehint!(values, nrow(stations_df))

        for i in 1:nrow(stations_df)
            x_idx = round(Int, x_coords[i]) + 1  # Convert 0-based to 1-based
            y_idx = round(Int, y_coords[i]) + 1

            # Bounds checking
            if x_idx < 1 || x_idx > GRID_SIZE || y_idx < 1 || y_idx > GRID_SIZE
                push!(values, NaN)
                continue
            end

            # NCDatasets exposes the NetCDF y/x dimensions in Julia indexing
            # order as x/y. Navigation audits against gauges confirm [x, y].
            val = precip[x_idx, y_idx]
            quality_value = quality === nothing ? 0 : quality[x_idx, y_idx]

            # Check for fill value or invalid
            if ismissing(val) || ismissing(quality_value) ||
               val == fill_val || val < 0 || val > 30 || Int(quality_value) >= 2
                push!(values, NaN)
            else
                push!(values, Float64(val))
            end
        end

        return values
    finally
        close(ds)
    end
end

"""
    aggregate_hourly(file_groups::Dict)

Aggregate 15min precipitation files to hourly sums.
file_groups: Dict mapping hour_start_datetime -> list of 4 NC file paths
Returns: DataFrame with time, station_id, fy4b_mm_h
"""
function aggregate_hourly(file_groups, stations_df)
    results = DataFrame(time=DateTime[], station_id=String[], fy4b_mm_h=Float64[])
    coord_cache = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()

    hours = sort(collect(keys(file_groups)))

    for (hour_index, hour_start) in enumerate(hours)
        files = file_groups[hour_start]

        starts = [first(parse_nc_filename(basename(file))) for file in files]
        minutes = sort(minute.(starts))
        @assert length(files) == 4 && minutes == [0, 15, 30, 45] "Non-complete hour reached aggregation: $hour_start"

        # The 15-minute QPE fields are rainfall rates. Their four-value mean is
        # numerically the one-hour accumulation in mm (sum(rate * 0.25 h)).
        hourly_sum = zeros(Float64, nrow(stations_df))
        has_missing = falses(nrow(stations_df))

        for nc_file in files
            vals = try
                extract_precipitation(nc_file, stations_df, coord_cache)
            catch error_value
                error("Failed to extract $(nc_file): $(sprint(showerror, error_value))")
            end
            for i in 1:length(vals)
                if isnan(vals[i])
                    has_missing[i] = true
                else
                    hourly_sum[i] += vals[i]
                end
            end
        end

        # Mark as NaN if any 15min value was missing
        for i in 1:length(hourly_sum)
            if has_missing[i]
                hourly_sum[i] = NaN
            else
                hourly_sum[i] /= 4.0
            end
        end

        # Add to results
        hour_end = aligned_hour_end(hour_start)
        time_str = Dates.format(hour_end, dateformat"yyyy-mm-ddTHH:MM:SS") * "Z"

        for (i, row) in enumerate(eachrow(stations_df))
            push!(results, (DateTime(replace(time_str, "Z" => "")), string(row.station_id), hourly_sum[i]))
        end
        if hour_index % 500 == 0 || hour_index == length(hours)
            println("  Aggregated $hour_index / $(length(hours)) strict hours")
        end
    end

    return results
end

# ---------------------------------------------------------------------------------------------
# Rebuild variants: a neighbourhood window and a parallax shift, sampled in one pass. None of the
# choices is tuned on gauges, so the rebuilt series can enter cross-validation without leaking.
# ---------------------------------------------------------------------------------------------

"""The 15-minute product's ceiling: it declares `valid_range = [0, 20]` mm/h."""
const FY4B_CAP = 20.0

"""
`true` when a raw value and its quality flag pass `extract_precipitation`'s rule. The `> 30`
bound is what drops the 65533 and 65535 sentinels: the files name their fill attribute
`FillValue`, not CF `_FillValue`, so nothing masks them on read, and valid values stop at 20.
"""
accepts_value(value, quality) =
    !ismissing(value) && !ismissing(quality) && value != Float32(65534.0) &&
    value >= 0 && value <= 30 && Int(quality) < 2

"""
Value of the `window`×`window` box centred on disk pixel `(ix, iy)` (1-based), as
`(value, cap_fraction)`; `block[1, 1]` is disk pixel `(x0, y0)` and `quality` its DQF block.

`window == 1` is the centre pixel exactly as `extract_precipitation` reads it. A wider window
averages the accepted pixels and needs more than half of them (5 of 9 for 3×3), a rule fixed in
advance rather than tuned; otherwise it is `NaN`. `cap_fraction` is the share of accepted pixels
at `FY4B_CAP`, where the product truncates heavier rain.
"""
function _window_value(block, quality, ix, iy, window; x0=1, y0=1)
    half = window ÷ 2
    total = 0.0
    accepted = 0
    capped = 0
    for jy in (iy - half):(iy + half), jx in (ix - half):(ix + half)
        bx = jx - x0 + 1
        by = jy - y0 + 1
        (1 <= bx <= size(block, 1) && 1 <= by <= size(block, 2)) || continue
        value = block[bx, by]
        accepts_value(value, quality === nothing ? 0 : quality[bx, by]) || continue
        total += Float64(value)
        accepted += 1
        Float64(value) >= FY4B_CAP && (capped += 1)
    end
    needed = window == 1 ? 1 : window^2 ÷ 2 + 1
    accepted >= needed || return (NaN, NaN)
    return (window == 1 ? total : total / accepted, capped / accepted)
end

_pixel(x, y) = isnan(x) || isnan(y) ? (0, 0) : (round(Int, x) + 1, round(Int, y) + 1)

"""
Every variant at every station from one FY4B segment, as
`(values, cap_fraction, sat_lon, software, shifted)` with the matrices `[station, variant]`.

`variants` holds `(name, window, height)`: the box size in pixels and the parallax cloud-top
height in metres (0 = none). The file is opened once and a single block covering every needed
pixel is read, so extra variants cost almost nothing. Grid coordinates are cached per
`(sat_lon, height)`, the satellite longitude coming from the file as in `extract_precipitation`.

`offsets = (; window, radius)` also returns `shifted[station, dx, dy]`: the `window` box value
with the station's pixel moved by every `(dx, dy)` in `-radius:radius` (x east-positive, y
south-positive, as image rows run). A later step picks a shift on training gauges only; this
function only tabulates the candidates. Without `offsets`, `shifted` is `nothing`.
"""
function sample_variants(nc_file, stations_df, variants, coord_cache; offsets=nothing)
    ds = Dataset(nc_file)
    try
        sat_lon = dataset_satellite_lon(ds, basename(nc_file))
        software = string(get(ds.attrib, "Version Of Software", ""))
        station_centres(height) = begin
            coords = get!(coord_cache, (sat_lon, Float64(height))) do
                [latlon_to_xy(row.lat, row.lon, sat_lon; height) for row in eachrow(stations_df)]
            end
            [_pixel(x, y) for (x, y) in coords]
        end
        centres = [station_centres(variant.height) for variant in variants]
        base = offsets === nothing ? nothing : station_centres(0.0)
        reach = maximum(variant.window ÷ 2 for variant in variants)
        offsets === nothing || (reach = max(reach, offsets.radius + offsets.window ÷ 2))
        all_centres = base === nothing ? centres : vcat(centres, [base])
        inside = [c for cs in all_centres for c in cs if c != (0, 0)]
        n, nv = nrow(stations_df), length(variants)
        values = fill(NaN, n, nv)
        caps = fill(NaN, n, nv)
        side = offsets === nothing ? 0 : 2 * offsets.radius + 1
        shifted = offsets === nothing ? nothing : fill(NaN, n, side, side)
        isempty(inside) && return values, caps, sat_lon, software, shifted
        x1 = clamp(minimum(first, inside) - reach, 1, GRID_SIZE)
        x2 = clamp(maximum(first, inside) + reach, 1, GRID_SIZE)
        y1 = clamp(minimum(last, inside) - reach, 1, GRID_SIZE)
        y2 = clamp(maximum(last, inside) + reach, 1, GRID_SIZE)
        # NCDatasets exposes the NetCDF y/x dimensions as [x, y], as in extract_precipitation.
        block = ds["Precipitation"][x1:x2, y1:y2]
        quality = haskey(ds, "DQF") ? ds["DQF"][x1:x2, y1:y2] : nothing
        for (k, variant) in enumerate(variants), i in 1:n
            ix, iy = centres[k][i]
            values[i, k], caps[i, k] =
                _window_value(block, quality, ix, iy, variant.window; x0=x1, y0=y1)
        end
        if offsets !== nothing
            r = offsets.radius
            for (b, dy) in enumerate(-r:r), (a, dx) in enumerate(-r:r), i in 1:n
                ix, iy = base[i]
                (ix, iy) == (0, 0) && continue
                shifted[i, a, b] = first(_window_value(
                    block, quality, ix + dx, iy + dy, offsets.window; x0=x1, y0=y1))
            end
        end
        return values, caps, sat_lon, software, shifted
    finally
        close(ds)
    end
end

"""
Hourly series for every variant from strict-complete hours, as
`(; times, values, cap_fraction, sat_lon, software)`.

`values` and `cap_fraction` map each variant name to a `[time, station]` matrix; `times` are the
hour-ending Beijing labels of `aggregate_hourly`, and `sat_lon`/`software` describe each hour's
first segment, from which the product regimes are read. The arithmetic is `aggregate_hourly`'s -
the four rates summed in file order and divided by 4, one missing segment voiding the
station-hour - so the `window == 1`, `height == 0` variant reproduces it bit for bit.

With `offsets` (see `sample_variants`), each hour's shifted grid `[station, dx, dy]`, aggregated
the same way, is handed to `on_shifted(hour_index, grid)` instead of being kept: at 169 shifts
the whole archive would not fit in memory.
"""
function aggregate_hourly_variants(
    file_groups, stations_df, variants;
    offsets=nothing, on_shifted=(_, _) -> nothing, progress_every::Int=500,
)
    hours = sort(collect(keys(file_groups)))
    n, nv = nrow(stations_df), length(variants)
    values = [fill(NaN, length(hours), n) for _ in 1:nv]
    caps = [fill(NaN, length(hours), n) for _ in 1:nv]
    sat_lons = fill(NaN, length(hours))
    software = fill("", length(hours))
    coord_cache = Dict{Tuple{Float64, Float64}, Vector{Tuple{Float64, Float64}}}()
    side = offsets === nothing ? 0 : 2 * offsets.radius + 1

    for (h, hour_start) in enumerate(hours)
        files = file_groups[hour_start]
        minutes = sort([minute(first(parse_nc_filename(basename(file)))) for file in files])
        @assert length(files) == 4 && minutes == [0, 15, 30, 45] "Non-complete hour reached aggregation: $hour_start"
        sums = zeros(n, nv)
        cap_sums = zeros(n, nv)
        gap = falses(n, nv)
        shifted_sums = zeros(n, side, side)
        shifted_gap = falses(n, side, side)
        for (segment, nc_file) in enumerate(files)
            segment_values, segment_caps, sat_lon, version, shifted = try
                sample_variants(nc_file, stations_df, variants, coord_cache; offsets)
            catch error_value
                error("Failed to extract $(nc_file): $(sprint(showerror, error_value))")
            end
            if shifted !== nothing
                for index in eachindex(shifted)
                    isnan(shifted[index]) ? (shifted_gap[index] = true) :
                        (shifted_sums[index] += shifted[index])
                end
            end
            if segment == 1
                sat_lons[h] = sat_lon
                software[h] = version
            end
            for k in 1:nv, i in 1:n
                if isnan(segment_values[i, k])
                    gap[i, k] = true
                else
                    sums[i, k] += segment_values[i, k]
                    cap_sums[i, k] += segment_caps[i, k]
                end
            end
        end
        for k in 1:nv, i in 1:n
            gap[i, k] && continue
            values[k][h, i] = sums[i, k] / 4.0
            caps[k][h, i] = cap_sums[i, k] / 4.0
        end
        if offsets !== nothing
            on_shifted(h, ifelse.(shifted_gap, NaN, shifted_sums ./ 4.0))
        end
        if progress_every > 0 && (h % progress_every == 0 || h == length(hours))
            println("  Sampled $h / $(length(hours)) strict hours")
        end
    end

    variant_names = [variant.name for variant in variants]
    return (;
        times=aligned_hour_end.(hours),
        values=Dict(zip(variant_names, values)),
        cap_fraction=Dict(zip(variant_names, caps)),
        sat_lon=sat_lons,
        software,
    )
end

"""Hours per time chunk of the shifted grid: a 31-day month, so one shift reads in ~36 chunks."""
const OFFSET_GRID_BLOCK = 744

"""
Writer for `fy4b_offset_grid.nc`, as `(on_shifted, finish, abort)`.

The variable `fy4b[offset_x, offset_y, station, time]` (Float32, deflated) is chunked one shift
by one month, so a later step can read any single shift's full series cheaply. Hours arrive one
at a time from `aggregate_hourly_variants` and are written a chunk-aligned block at a time, so
each chunk is compressed once. The file appears under its name only after `finish()`.
"""
function _offset_grid_writer(path, offsets, station_ids, times)
    r = offsets.radius
    side = 2r + 1
    n, total = length(station_ids), length(times)
    temporary = string(path, ".tmp-", getpid())
    mkpath(dirname(path))
    ds = NCDataset(temporary, "c")
    defDim(ds, "offset_x", side)
    defDim(ds, "offset_y", side)
    defDim(ds, "station", n)
    defDim(ds, "time", total)
    defVar(ds, "offset_x", collect(-r:r), ("offset_x",);
        attrib=["units" => "pixels", "long_name" => "sample pixel moved east (+) or west (-)"])
    defVar(ds, "offset_y", collect(-r:r), ("offset_y",);
        attrib=["units" => "pixels",
            "long_name" => "sample pixel moved south (+) or north (-), as image rows run"])
    defVar(ds, "station_id", station_ids, ("station",))
    defVar(ds, "time_label", Dates.format.(times, dateformat"yyyy-mm-ddTHH:MM:SS") .* "Z",
        ("time",); attrib=["long_name" => "hour-ending Beijing label, as in the wide CSVs"])
    grid = defVar(ds, "fy4b", Float32, ("offset_x", "offset_y", "station", "time");
        deflatelevel=4, shuffle=true, chunksizes=[1, 1, n, min(OFFSET_GRID_BLOCK, total)],
        attrib=["units" => "mm/h", "window" => offsets.window, "radius" => r,
            "long_name" => "hourly FY4B, window mean around the shifted station pixel",
            # No fill value: a gap is stored and read back as NaN, as in the wide CSVs.
            "missing_data" => "NaN"])
    buffer = fill(NaN32, side, side, n, OFFSET_GRID_BLOCK)
    filled = 0
    start = 1
    flush!() = begin
        if filled > 0
            grid[:, :, :, start:(start + filled - 1)] = buffer[:, :, :, 1:filled]
            start += filled
            filled = 0
        end
    end
    on_shifted(hour_index, hour_grid) = begin
        hour_index == start + filled || error("shifted grids must arrive in hour order")
        buffer[:, :, :, filled + 1] = permutedims(hour_grid, (2, 3, 1))
        filled += 1
        filled == OFFSET_GRID_BLOCK && flush!()
    end
    finish() = begin
        flush!()
        close(ds)
        mv(temporary, path; force=true)
    end
    abort() = begin
        close(ds)
        rm(temporary; force=true)
    end
    return on_shifted, finish, abort
end

"""Wide `time, <station>...` table in the layout `aggregate_fy4b_hourly` writes."""
function _wide_table(time_text, station_ids, matrix)
    table = DataFrame(time=time_text)
    for (column, station_id) in enumerate(station_ids)
        table[!, station_id] = matrix[:, column]
    end
    return table
end

"""
    aggregate_fy4b_variants(; kwargs...)

Rebuild the hourly FY4B station series in several variants from one pass over the raw tree.

Hour selection is `aggregate_fy4b_hourly`'s: strict-complete hours only, from 15-minute files.
Writes to `outdir`:
- `fy4b_<variant>.csv`: the wide series per variant;
- `fy4b_<variant>_capfrac.csv`: per station-hour, the share of accepted pixels at the 20 mm/h
  cap, so a later step can treat them as censored;
- `fy4b_hourly_qc.csv`: the per-hour completeness audit;
- `fy4b_hour_meta.csv`: each hour's satellite longitude and software version;
- with `offset_grid = (; window, radius)`, `fy4b_offset_grid.nc`: the `window` box value at
  every pixel shift up to `radius`, for a shift chosen later on training gauges only.
"""
function aggregate_fy4b_variants(;
    data_dir::AbstractString=FY4B_DATA_DIR,
    station_meta_file::AbstractString,
    outdir::AbstractString,
    variants,
    years=2022:2024,
    months=1:12,
    offset_grid=nothing,
)
    variant_names = [variant.name for variant in variants]
    allunique(variant_names) || error("Variant names must be unique: $variant_names")
    all(variant -> isodd(variant.window) && variant.window >= 1 && variant.height >= 0, variants) ||
        error("Variant windows must be odd and positive, heights non-negative")
    offset_grid === nothing || (isodd(offset_grid.window) && offset_grid.radius >= 0) ||
        error("The offset grid needs an odd window and a non-negative radius")

    stations_df = load_stations(station_meta_file)
    all(name -> name in names(stations_df), ["station_id", "lon", "lat"]) ||
        error("Station metadata must contain station_id, lon, and lat columns")
    nrow(unique(stations_df, :station_id)) == nrow(stations_df) ||
        error("Station metadata contains duplicate station_id values")
    println("Loaded $(nrow(stations_df)) stations; variants: $(join(variant_names, ", "))")

    selected_months = normalize_month_selection(months, 1, 13)
    file_groups = find_nc_files(data_dir; years, months=selected_months)
    unreadable_files = find_unreadable_files(file_groups)
    for (file, message) in sort(collect(unreadable_files); by=first)
        @warn "Unreadable FY4B segment; treating as missing" file message
    end
    native_hours = find_native_hourly_times(data_dir)
    qc_df, complete_groups = build_hourly_qc(
        file_groups, native_hours; unreadable_files, years, months=selected_months,
    )
    write_csv_atomic(joinpath(outdir, "fy4b_hourly_qc.csv"), qc_df)
    complete_count = length(complete_groups)
    println("Strict complete hours: $complete_count of $(nrow(qc_df))")
    complete_count > 0 || error("No strict complete FY4B hours were found; QC is in $outdir")

    station_ids = string.(stations_df.station_id)
    result = if offset_grid === nothing
        aggregate_hourly_variants(complete_groups, stations_df, variants)
    else
        hours = sort(collect(keys(complete_groups)))
        on_shifted, finish, abort = _offset_grid_writer(
            joinpath(outdir, "fy4b_offset_grid.nc"), offset_grid, station_ids,
            aligned_hour_end.(hours),
        )
        try
            sampled = aggregate_hourly_variants(
                complete_groups, stations_df, variants; offsets=offset_grid, on_shifted,
            )
            finish()
            sampled
        catch
            abort()
            rethrow()
        end
    end
    time_text = Dates.format.(result.times, dateformat"yyyy-mm-ddTHH:MM:SS") .* "Z"
    for name in variant_names
        write_csv_atomic(joinpath(outdir, "fy4b_$(name).csv"),
            _wide_table(time_text, station_ids, result.values[name]))
        write_csv_atomic(joinpath(outdir, "fy4b_$(name)_capfrac.csv"),
            _wide_table(time_text, station_ids, result.cap_fraction[name]))
    end
    write_csv_atomic(joinpath(outdir, "fy4b_hour_meta.csv"),
        DataFrame(time=time_text, sat_lon=result.sat_lon, software=result.software))
    println("Wrote $(length(variant_names)) variants x $complete_count hours to $outdir")
    return (; outdir, complete_count, unreadable_files)
end

"""
`true` for a macOS AppleDouble companion (`._<name>`), which macOS writes beside each file
copied to a non-Apple volume such as the exFAT drive holding the raw tree. It carries the real
file's name, so it would pass the name filter - and sort ahead of the real file.
"""
is_appledouble(filename) = startswith(filename, "._")

"""
    find_nc_files(data_dir)

Find all NC files and group 15min files by hour.
Returns: Dict mapping hour_start -> list of NC file paths
"""
function find_nc_files(data_dir; years=2022:2025, months=6:9)
    file_groups = Dict{DateTime, Vector{String}}()
    candidates = Dict{Tuple{DateTime, Float64}, Vector{String}}()

    # Downloads are stored below per-year batch directories, not year/month.
    for (dir, _, files) in walkdir(data_dir), file in files
        is_appledouble(file) && continue
        endswith(uppercase(file), ".NC") || continue
        is_15min_file(file) || continue
        start_dt, _, sat_lon = parse_nc_filename(file)
        year(start_dt) in years || continue
        month(start_dt) in months || continue
        push!(get!(candidates, (start_dt, sat_lon), String[]), joinpath(dir, file))
    end

    duplicate_file_count = 0
    for ((start_dt, _), paths) in candidates
        sort!(paths)
        duplicate_file_count += length(paths) - 1
        hour_start = DateTime(year(start_dt), month(start_dt), day(start_dt), hour(start_dt))
        # Of several copies keep the first that opens, so a truncated download sorting first
        # cannot shadow a good copy elsewhere in the tree. A lone copy is kept as is.
        index = length(paths) == 1 ? 1 :
            something(findfirst(path -> segment_problem(path) === nothing, paths), 1)
        push!(get!(file_groups, hour_start, String[]), paths[index])
    end

    # Sort files within each hour group
    for (hour, files) in file_groups
        sort!(files)
    end
    duplicate_file_count > 0 && @info "Deduplicated repeated FY4B 15-minute files" duplicate_file_count

    return file_groups
end

"""`nothing` when `file` opens as a full-disk FY4B segment, otherwise the reason it does not."""
function segment_problem(file)
    try
        ds = Dataset(file)
        try
            haskey(ds, "Precipitation") || error("missing Precipitation variable")
            precipitation_size = size(ds["Precipitation"])
            precipitation_size == (GRID_SIZE, GRID_SIZE) ||
                error("unexpected Precipitation size $precipitation_size")
        finally
            close(ds)
        end
    catch error_value
        return sprint(showerror, error_value)
    end
    return nothing
end

"""Open every selected NetCDF and return the paths that fail structural validation."""
function find_unreadable_files(file_groups)
    files = sort(unique(vcat(values(file_groups)...)))
    unreadable = Dict{String, String}()
    for (index, file) in enumerate(files)
        problem = segment_problem(file)
        problem === nothing || (unreadable[file] = problem)
        if index % 5_000 == 0 || index == length(files)
            println("  Validated $index / $(length(files)) NetCDF files")
        end
    end
    return unreadable
end

"""Return hour starts for exact native one-hour FY4B products in the raw download tree."""
function find_native_hourly_times(data_dir)
    native_hours = Set{DateTime}()
    for (dir, _, files) in walkdir(data_dir), file in files
        is_appledouble(file) && continue
        endswith(file, ".NC") || continue
        dt_info = parse_nc_filename(file)
        dt_info === nothing && continue
        start_dt, end_dt, _ = dt_info
        if minute(start_dt) == 0 && second(start_dt) == 0 &&
           end_dt == start_dt + Hour(1) - Second(1)
            push!(native_hours, start_dt)
        end
    end
    return native_hours
end

function normalize_month_selection(months, start_month::Int, stop_month::Int)
    selected = months === nothing ?
        collect(start_month:(stop_month - 1)) :
        sort(unique(Int.(collect(months))))
    isempty(selected) && error("At least one month must be selected")
    all(month_value -> 1 <= month_value <= 12, selected) ||
        error("Months must be between 1 and 12: $selected")
    return selected
end

"""Build one auditable row for every expected hour in the selected year-months."""
function build_hourly_qc(
    file_groups, native_hours;
    unreadable_files=Dict{String, String}(), years=2022:2025,
    start_month::Int=6, stop_month::Int=10, months=nothing,
)
    rows = NamedTuple[]
    complete_groups = Dict{DateTime, Vector{String}}()
    expected_minutes = [0, 15, 30, 45]
    selected_months = normalize_month_selection(months, start_month, stop_month)

    for year_value in years
        for month_value in selected_months
            hour_start = DateTime(year_value, month_value, 1)
            stop = hour_start + Month(1)
            while hour_start < stop
                files = sort(get(file_groups, hour_start, String[]))
                readable_files = [file for file in files if !haskey(unreadable_files, file)]
                unreadable_hour_files = [file for file in files if haskey(unreadable_files, file)]
                starts = DateTime[]
                all_starts = DateTime[]
                for file in readable_files
                    info = parse_nc_filename(basename(file))
                    info === nothing || push!(starts, info[1])
                end
                for file in files
                    info = parse_nc_filename(basename(file))
                    info === nothing || push!(all_starts, info[1])
                end
                present_minutes = sort(unique(minute.(starts)))
                missing_minutes = setdiff(expected_minutes, present_minutes)
                duplicate_count = length(all_starts) - length(unique(minute.(all_starts)))
                complete = length(readable_files) == 4 &&
                    present_minutes == expected_minutes && duplicate_count == 0
                if complete
                    complete_groups[hour_start] = readable_files
                end

                present_text = join([lpad(string(value), 2, '0') for value in present_minutes], ";")
                missing_text = join([lpad(string(value), 2, '0') for value in missing_minutes], ";")
                source_text = join(basename.(files), ";")
                push!(rows, (
                    hour_start = Dates.format(hour_start, dateformat"yyyy-mm-ddTHH:MM:SS") * "Z",
                    hour_end = Dates.format(hour_start + Hour(1), dateformat"yyyy-mm-ddTHH:MM:SS") * "Z",
                    status = complete ? "complete" : "incomplete",
                    file_count = length(files),
                    segment_count = length(readable_files),
                    present_minutes = present_text,
                    missing_minutes = missing_text,
                    missing_count = length(missing_minutes),
                    duplicate_count = duplicate_count,
                    unreadable_count = length(unreadable_hour_files),
                    unreadable_files = join(basename.(unreadable_hour_files), ";"),
                    native_hourly_available = hour_start in native_hours,
                    source_files = source_text,
                ))
                hour_start += Hour(1)
            end
        end
    end

    return DataFrame(rows), complete_groups
end

"""
    pivot_to_wide(long_df)

Convert long format (time, station_id, fy4b_mm_h) to wide format
(time, station_1, station_2, ...)
"""
function pivot_to_wide(long_df)
    # Ensure station_id is string
    long_df.station_id = string.(long_df.station_id)

    # Pivot to wide format
    wide = unstack(long_df, :time, :station_id, :fy4b_mm_h)

    # Sort by time
    sort!(wide, :time)

    return wide
end

"""
    aggregate_fy4b_hourly(; kwargs...)

Main pipeline to aggregate FY4B data.
"""
function aggregate_fy4b_hourly(;
    data_dir::AbstractString=FY4B_DATA_DIR,
    station_meta_file::AbstractString=STATION_META_FILE,
    output_file::AbstractString=OUTPUT_FILE,
    qc_file::AbstractString=QC_FILE,
    years=2022:2025,
    start_month::Int=6,
    stop_month::Int=10,
    months=nothing,
)
    println("="^60)
    println("FY4B 15min -> Hourly Aggregation Pipeline")
    println("="^60)

    # Step 1: Load station metadata
    println("\n[1/5] Loading station metadata...")
    stations_df = load_stations(station_meta_file)
    all(name -> name in names(stations_df), ["station_id", "lon", "lat"]) ||
        error("Station metadata must contain station_id, lon, and lat columns")
    nrow(unique(stations_df, :station_id)) == nrow(stations_df) ||
        error("Station metadata contains duplicate station_id values")
    println("  Loaded $(nrow(stations_df)) stations")

    # Step 2: Find and group NC files
    println("\n[2/5] Finding 15min NC files...")
    selected_months = normalize_month_selection(months, start_month, stop_month)
    file_groups = find_nc_files(data_dir; years, months=selected_months)
    println("  Validating selected NetCDF structures...")
    unreadable_files = find_unreadable_files(file_groups)
    for (file, message) in sort(collect(unreadable_files); by=first)
        @warn "Unreadable FY4B segment; treating as missing" file message
    end
    native_hours = find_native_hourly_times(data_dir)
    qc_df, complete_groups = build_hourly_qc(
        file_groups, native_hours;
        unreadable_files, years, start_month, stop_month, months=selected_months,
    )
    total_hours = length(file_groups)
    total_files = sum(length.(values(file_groups)))
    complete_count = count(==("complete"), qc_df.status)
    incomplete_count = count(==("incomplete"), qc_df.status)
    missing_segment_count = sum(qc_df.missing_count)
    native_fallback_count = count((qc_df.status .== "incomplete") .& qc_df.native_hourly_available)
    println("  Found $total_hours represented hours with $total_files total 15min files")
    println("  Strict complete hours: $complete_count; incomplete hours: $incomplete_count")
    println("  Missing 15min segments: $missing_segment_count")
    println("  Incomplete hours with native hourly product (QC only): $native_fallback_count")
    @assert complete_count + incomplete_count == nrow(qc_df)
    @assert complete_count == length(complete_groups)
    write_csv_atomic(qc_file, qc_df)
    println("  QC saved: $qc_file")
    sat_lons = sort(unique(parse_satellite_lon(basename(file)) for files in values(file_groups) for file in files))
    sat_lon_text = join(string.(sat_lons), ", ")
    println("  Satellite subpoints in selected files: $sat_lon_text E")
    complete_count > 0 || error(
        "No strict complete FY4B hours were found; QC was written to $qc_file",
    )

    # Step 3: Aggregate to hourly
    println("\n[3/5] Aggregating precipitation to hourly...")
    long_df = aggregate_hourly(complete_groups, stations_df)
    println("  Aggregated $(nrow(long_df)) station-hour values")

    # Step 4: Convert to wide format
    println("\n[4/5] Converting to wide format...")
    wide_df = pivot_to_wide(long_df)
    println("  Wide format: $(size(wide_df, 1)) time rows × $(size(wide_df, 2)-1) station columns")

    # Step 5: Write output
    println("\n[5/5] Writing output...")
    # Format time column as ISO with Z
    wide_df.time = Dates.format.(wide_df.time, dateformat"yyyy-mm-ddTHH:MM:SS") .* "Z"
    write_csv_atomic(output_file, wide_df)
    println("  Saved: $output_file")

    # Summary statistics
    println("\n" * "="^60)
    println("Summary:")
    println("  Output dimensions: $(size(wide_df))")
    println("  Time range: $(wide_df.time[1]) to $(wide_df.time[end])")
    println("  Strict complete rows: $complete_count")
    println("  Station columns: $(size(wide_df, 2)-1)")
    println("="^60)
    return (; output_file, qc_file, complete_count, incomplete_count, unreadable_files)
end

end # module
