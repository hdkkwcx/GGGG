"""
FY4B's Cloud Type product (`CLT`) at the gauges, for screening a false-rain filter.

`CLT` is the AGRI L2 cloud type: full disk, 4 km, every 15 minutes. It is on the same nominal grid
and scans as the QPE files, and NCDatasets exposes both as `[x, y]`, so a gauge's QPE pixel is its
CLT pixel. The class codes are `CLT_CLASSES`; 126 (space) and 127 (fill) are not classes. DQF is
not read, because class 0 already is the cloud mask's clear call.

- `find_clt_files`: the 15-minute CLT files under a tree, grouped by UTC hour, with a QC table.
- `clt_class_counts`: per station, pixel shift and complete hour, how often each class occurs
  over the `window` box and the hour's four quarters.
- `dominant_class`: each cell's most frequent class.

Hours follow the FY4B rebuild (`FY4BPreprocessing`). They are strict, so all four quarters or no
value, and labelled by the end of the hour in Beijing time (`aligned_hour_end`).
"""
module FY4BCloudType

using DataFrames, Dates, NCDatasets
using Main.FY4BPreprocessing
using Main: FY4BPreprocessing

export CLT_CLASSES, CLT_CLASS_NAMES, find_clt_files, clt_class_counts, dominant_class

"""The `CLT` codes that are cloud-type classes, in the order of `clt_class_counts`' last axis."""
const CLT_CLASSES = (0, 2, 3, 4, 5, 6, 7, 9)
const CLT_CLASS_NAMES =
    ("clear", "water", "supercooled", "mixed", "ice", "cirrus", "overlap", "uncertain")
const QUARTER_MINUTES = (0, 15, 30, 45)

"""Index of a raw `CLT` code in `CLT_CLASSES`, or 0 for space, fill and anything else."""
_class_index(code::Integer) = something(findfirst(==(code), CLT_CLASSES), 0)
const CLASS_INDEX = Tuple(_class_index(code) for code in 0:127)
_class_of(code) = 0 <= code <= 127 ? CLASS_INDEX[code + 1] : 0

"""`true` for a 15-minute CLT file name (not a QPE file, not a macOS `._` companion)."""
is_clt_file(filename) = occursin("_L2-_CLT-_", filename) &&
    endswith(uppercase(filename), ".NC") && !FY4BPreprocessing.is_appledouble(filename) &&
    FY4BPreprocessing.is_15min_file(filename)

"""`nothing` when `file` opens with a full-disk `CLT` variable, otherwise the reason it does not."""
function clt_problem(file)
    try
        NCDataset(file) do ds
            haskey(ds, "CLT") || return "no CLT variable"
            size(ds["CLT"]) == (FY4BPreprocessing.GRID_SIZE, FY4BPreprocessing.GRID_SIZE) ||
                return "CLT is $(size(ds["CLT"]))"
            return nothing
        end
    catch err
        return sprint(showerror, err)
    end
end

"""
The 15-minute CLT files under `dir`, as `(groups, qc)`.

- `groups`: UTC hour start => the hour's files, sorted. Copies of one quarter (same start and
  satellite position) are deduplicated: the first that opens is kept, as `find_nc_files` does.
- `qc`: one row per UTC hour of `years` x `months`, with its Beijing label `time`, the quarters
  present and missing, and `complete` (all four quarters, one file each).
"""
function find_clt_files(dir; years, months)
    candidates = Dict{Tuple{DateTime,Float64},Vector{String}}()
    for (root, _, files) in walkdir(dir), file in files
        is_clt_file(file) || continue
        start, _, sat_lon = parse_nc_filename(file)
        year(start) in years && month(start) in months || continue
        push!(get!(candidates, (start, sat_lon), String[]), joinpath(root, file))
    end
    groups = Dict{DateTime,Vector{String}}()
    for ((start, _), paths) in candidates
        sort!(paths)
        index = length(paths) == 1 ? 1 :
            something(findfirst(path -> clt_problem(path) === nothing, paths), 1)
        push!(get!(groups, floor(start, Hour), String[]), paths[index])
    end
    foreach(sort!, values(groups))

    rows = NamedTuple[]
    for y in years, m in months
        first_hour = DateTime(y, m, 1)
        for hour_start in first_hour:Hour(1):(first_hour + Month(1) - Hour(1))
            files = get(groups, hour_start, String[])
            minutes = [minute(first(parse_nc_filename(basename(file)))) for file in files]
            present = sort(unique(minutes))
            missing_minutes = [q for q in QUARTER_MINUTES if !(q in present)]
            push!(rows, (; hour_start, time=FY4BPreprocessing.aligned_hour_end(hour_start),
                file_count=length(files), present_minutes=join(present, ";"),
                missing_minutes=join(missing_minutes, ";"),
                complete=_complete(files), files=join(basename.(files), ";")))
        end
    end
    return groups, DataFrame(rows)
end

"""All four quarters, each from exactly one file."""
_complete(files) = length(files) == 4 &&
    sort([minute(first(parse_nc_filename(basename(file)))) for file in files]) ==
    collect(QUARTER_MINUTES)

"""
Class counts at every station, as `(; counts, hour_start, time)`.

- `counts[station, s, hour, c]` is how many of the hour's pixel-quarters in the `window` x `window`
  box hold class `CLT_CLASSES[c]`. The box is centred on the station's nominal pixel moved by
  `shifts[s] = (dx, dy)`, x east-positive and y south-positive, as in `sample_variants`' offset
  grid. Pixels holding space, fill or an unknown code, or lying off the disk, are not counted.
- Only complete hours (all four quarters) are kept: `hour_start` (UTC) and `time` (the Beijing
  hour-end label) list them in order.

`stations` needs `lat` and `lon`. A station's pixel follows each file's satellite position, as
the QPE rebuild does. Each file is opened once, and one block around every station is read.
"""
function clt_class_counts(groups, stations, shifts::AbstractVector{<:Tuple{Int,Int}};
                          window::Int=3)
    hours = sort([hour for (hour, files) in groups if _complete(files)])
    counts = zeros(UInt8, nrow(stations), length(shifts), length(hours), length(CLT_CLASSES))
    centres = Dict{Float64,Vector{Tuple{Int,Int}}}()
    reach = window ÷ 2 + maximum(s -> max(abs(s[1]), abs(s[2])), shifts)
    size_ = FY4BPreprocessing.GRID_SIZE
    for (h, hour) in enumerate(hours), file in groups[hour]
        NCDataset(file) do ds
            sat_lon = FY4BPreprocessing.dataset_satellite_lon(ds, basename(file))
            pixels = get!(centres, sat_lon) do
                [FY4BPreprocessing._pixel(FY4BPreprocessing.latlon_to_xy(row.lat, row.lon,
                    sat_lon)...) for row in eachrow(stations)]
            end
            inside = filter(!=((0, 0)), pixels)
            isempty(inside) && return
            x1 = clamp(minimum(first, inside) - reach, 1, size_)
            x2 = clamp(maximum(first, inside) + reach, 1, size_)
            y1 = clamp(minimum(last, inside) - reach, 1, size_)
            y2 = clamp(maximum(last, inside) + reach, 1, size_)
            # The raw Int8 codes: the files' `FillValue` is not CF `_FillValue`, so nothing is
            # masked either way, and `.var` skips `_Unsigned` decoding.
            block = ds["CLT"].var[x1:x2, y1:y2]
            half = window ÷ 2
            @inbounds for i in eachindex(pixels), (s, (dx, dy)) in enumerate(shifts)
                pixels[i] == (0, 0) && continue
                cx, cy = pixels[i][1] + dx, pixels[i][2] + dy
                for x in (cx - half):(cx + half), y in (cy - half):(cy + half)
                    x1 <= x <= x2 && y1 <= y <= y2 || continue
                    c = _class_of(Int(block[x - x1 + 1, y - y1 + 1]))
                    c == 0 || (counts[i, s, h, c] += 0x01)
                end
            end
        end
    end
    return (; counts, hour_start=hours, time=FY4BPreprocessing.aligned_hour_end.(hours))
end

"""
The most frequent class of every `[station, shift, hour]` cell, as a `CLT_CLASSES` code, or -1
when fewer than `min_valid` pixel-quarters hold a class. A tie goes to the smaller code.
"""
function dominant_class(counts::AbstractArray{<:Integer,4}; min_valid::Int)
    out = fill(Int8(-1), size(counts, 1), size(counts, 2), size(counts, 3))
    @inbounds for h in axes(counts, 3), s in axes(counts, 2), i in axes(counts, 1)
        total = 0
        best, best_count = 0, -1
        for c in axes(counts, 4)
            n = Int(counts[i, s, h, c])
            total += n
            n > best_count && ((best, best_count) = (c, n))
        end
        total >= min_valid && (out[i, s, h] = Int8(CLT_CLASSES[best]))
    end
    return out
end

end # module
