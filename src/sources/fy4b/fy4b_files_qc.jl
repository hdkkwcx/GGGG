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
