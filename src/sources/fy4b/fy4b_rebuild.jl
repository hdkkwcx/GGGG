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
