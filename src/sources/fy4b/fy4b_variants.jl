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
