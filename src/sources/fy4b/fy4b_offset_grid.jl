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
