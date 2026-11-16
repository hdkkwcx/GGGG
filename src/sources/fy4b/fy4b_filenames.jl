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
