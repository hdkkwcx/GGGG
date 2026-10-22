using Test, DataFrames, Dates, NCDatasets

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("IMERGUncalStations")
const IUS = Main.IMERGUncalStations

@testset "imerg_cell picks the cells the original IMERG Uncal series used" begin
    # Gauges on a cell edge, with the cell each one's published series matches (global 0-based
    # index; centre = first centre + 0.1 * index). Ties do not all fall the same way.
    @test IUS.imerg_cell(110.6, -179.95) == 2905     # centre 110.55, west of the edge
    @test IUS.imerg_cell(110.0, -179.95) == 2899     # 109.95
    @test IUS.imerg_cell(110.5, -179.95) == 2904     # 110.45
    @test IUS.imerg_cell(32.4, -89.95) == 1223       # 32.35, south of the edge
    @test IUS.imerg_cell(32.3, -89.95) == 1223       # 32.35, north of the edge
    @test IUS.imerg_cell(32.8, -89.95) == 1228       # 32.85
    @test IUS.imerg_cell(110.233, -179.95) == 2902   # an interior gauge: its nearest centre, 110.25
    # The first centre of the grid the GES DISC subsets start at (lon index 2886, lat index 1192).
    @test IUS.imerg_cell(108.649994f0, -179.95) == 2886
    @test IUS.imerg_cell(29.25f0, -89.95) == 1192
end

"""Write one synthetic half-hourly subset: `rate[lat, lon]` in mm/h, `missing` for a fill value."""
function write_halfhour(dir, start::DateTime, rate)
    daydir = joinpath(dir, string(year(start)), lpad(dayofyear(start), 3, '0'))
    mkpath(daydir)
    stop = start + Minute(29) + Second(59)
    name = "3B-HHR.MS.MRG.3IMERG.$(Dates.format(start, "yyyymmdd"))-S$(Dates.format(start, "HHMMSS"))-" *
        "E$(Dates.format(stop, "HHMMSS")).0000.V07B.subset.nc4"
    NCDataset(joinpath(daydir, name), "c") do ds
        defVar(ds, "lon", Float32[110.05, 110.15, 110.25, 110.35], ("lon",))
        defVar(ds, "lat", Float32[31.05, 31.15, 31.25], ("lat",))
        defVar(ds, "precipitationUncal", reshape(rate, 3, 4, 1), ("lat", "lon", "time");
            fillvalue=-9999.9f0)
    end
end

@testset "imerg_hourly_station_series averages the two half hours at each gauge's cell" begin
    dir = mktempdir()
    # A rate that identifies its half hour and its cell, so every hourly value is predictable.
    rate(t, i, j) = Float32(hour(t) + minute(t) / 60 + 0.1i + 0.01j)
    days = Date(2023, 2, 1):Day(1):Date(2023, 2, 28)
    for day in days, k in 0:47
        t = DateTime(day) + Minute(30k)
        grid = Union{Missing,Float32}[rate(t, i, j) for i in 1:3, j in 1:4]
        t == DateTime(2023, 2, 1, 0, 30) && (grid[2, 2] = missing)
        write_halfhour(dir, t, grid)
    end

    files = IUS.imerg_halfhour_files(dir, Date(2023, 2, 3))
    @test length(files) == 48
    @test first(files)[1] == DateTime(2023, 2, 3) && last(files)[1] == DateTime(2023, 2, 3, 23, 30)
    @test issorted(first.(files))

    stations = DataFrame(station_id=["a", "b", "far"], lon=[110.15, 110.33, 112.0], lat=[31.15, 31.06, 31.2])
    series = IUS.imerg_hourly_station_series(dir, stations, 2023, 2)
    @test series.outside == ["far"]
    table = series.table
    @test names(table) == ["time", "station_id", "gpm_mm_h"]
    @test nrow(table) == 28 * 24 * 2
    value(id, time) = only(table.gpm_mm_h[(table.station_id .== id) .& (table.time .== time)])
    t = DateTime(2023, 2, 14, 5)
    # Station a reads cell (2, 2), station b cell (1, 4): the mean of the two half-hour rates,
    # each widened to Float64 first.
    @test value("a", "2023-02-14T05:00:00Z") == (Float64(rate(t, 2, 2)) + Float64(rate(t + Minute(30), 2, 2))) / 2
    @test value("b", "2023-02-14T05:00:00Z") == (Float64(rate(t, 1, 4)) + Float64(rate(t + Minute(30), 1, 4))) / 2
    # A filled half hour makes its hour NaN, at that cell only.
    @test isnan(value("a", "2023-02-01T00:00:00Z"))
    @test !isnan(value("b", "2023-02-01T00:00:00Z"))

    # A day with a missing file is an error, not a silently short day.
    rm(last(IUS.imerg_halfhour_files(dir, Date(2023, 2, 10)))[2])
    @test_throws ErrorException IUS.imerg_halfhour_files(dir, Date(2023, 2, 10))
end
