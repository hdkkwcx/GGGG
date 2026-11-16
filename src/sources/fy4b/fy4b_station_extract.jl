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
