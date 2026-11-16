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
