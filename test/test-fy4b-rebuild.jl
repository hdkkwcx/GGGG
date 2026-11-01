using Test, DataFrames, Dates, NCDatasets

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("FY4BPreprocessing")
const FY4B = Main.FY4BPreprocessing

# A point in the middle of the study box, seen from both FY4B positions.
const STUDY_LON, STUDY_LAT = 110.5, 32.3

"""Great-circle distance (km) and initial bearing (degrees from north) from `a` to `b`."""
function distance_bearing(a, b)
    (lat1, lon1), (lat2, lon2) = deg2rad.(a), deg2rad.(b)
    h = sin((lat2 - lat1) / 2)^2 + cos(lat1) * cos(lat2) * sin((lon2 - lon1) / 2)^2
    distance = 2 * 6371.0 * asin(sqrt(h))
    bearing = rad2deg(atan(sin(lon2 - lon1) * cos(lat2),
        cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon2 - lon1)))
    return distance, mod(bearing, 360)
end

"""The surface point whose grid position is `(x, y)`, by a fine search around the study point."""
function ground_point(x, y, sat_lon)
    best, where = Inf, (NaN, NaN)
    for lat in (STUDY_LAT - 0.2):0.001:(STUDY_LAT + 0.2), lon in (STUDY_LON - 0.2):0.001:(STUDY_LON + 0.2)
        gx, gy = FY4B.latlon_to_xy(lat, lon, sat_lon)
        miss = (gx - x)^2 + (gy - y)^2
        miss < best && ((best, where) = (miss, (lat, lon)))
    end
    return where
end

@testset "AppleDouble companions never stand in for a segment" begin
    mktempdir() do dir
        name = "FY4B-_AGRI--_N_DISK_1330E_L2-_QPE-_MULT_NOM_20220601000000_20220601001459_4000M_V0001.NC"
        touch(joinpath(dir, name))
        touch(joinpath(dir, "._" * name))          # sorts first: it would win the dedup
        groups = FY4B.find_nc_files(dir; years=2022:2022, months=6:6)
        @test only(groups[DateTime(2022, 6, 1, 0)]) == joinpath(dir, name)
        @test FY4B.is_appledouble("._" * name) && !FY4B.is_appledouble(name)
    end
end

@testset "parallax lifts the sample point along the line of sight" begin
    for sat_lon in (133.0, 105.0)
        @test FY4B.latlon_to_xy(STUDY_LAT, STUDY_LON, sat_lon; height=0.0) ==
            FY4B.latlon_to_xy(STUDY_LAT, STUDY_LON, sat_lon)
    end
    # A 10 km cloud top over the study point shows up where this surface point does: ~9.9 km to
    # the north-west seen from 133E (zenith ~45 deg), ~7.8 km to the north from 105E (~38 deg).
    for (sat_lon, expected_km, bearing_range) in ((133.0, 9.9, (290, 340)), (105.0, 7.8, (0, 25)))
        x, y = FY4B.latlon_to_xy(STUDY_LAT, STUDY_LON, sat_lon; height=10_000.0)
        distance, bearing = distance_bearing((STUDY_LAT, STUDY_LON), ground_point(x, y, sat_lon))
        @test isapprox(distance, expected_km; atol=0.3)
        @test bearing_range[1] <= bearing <= bearing_range[2]
    end
end

@testset "window sampling keeps extract_precipitation's QC" begin
    block = Float32[1 2 3; 4 5 6; 7 8 20]
    good = zeros(Int8, 3, 3)
    @test FY4B._window_value(block, good, 2, 2, 1) == (5.0, 0.0)          # the centre, as read
    value, cap = FY4B._window_value(block, good, 2, 2, 3)
    @test value ≈ sum(Float64.(block)) / 9 && cap ≈ 1 / 9                  # one pixel at the cap
    # Rejected: DQF >= 2, the 65534 fill, and anything above 30 (the 65533/65535 sentinels).
    rejected = copy(block)
    rejected[1, 1] = 65534; rejected[1, 2] = 65535
    flags = copy(good)
    flags[1, 3] = 2
    value, _ = FY4B._window_value(rejected, flags, 2, 2, 3)
    @test value ≈ (4 + 5 + 6 + 7 + 8 + 20) / 6
    flags[2, 1] = 3                                                         # 5 accepted: still enough
    @test !isnan(first(FY4B._window_value(rejected, flags, 2, 2, 3)))
    flags[2, 2] = 2                                                         # 4 of 9: not enough
    @test all(isnan, FY4B._window_value(rejected, flags, 2, 2, 3))
    @test all(isnan, FY4B._window_value(rejected, flags, 2, 2, 1))          # rejected centre
    # A block offset to disk pixel (10, 20), and a box running off its edge.
    @test FY4B._window_value(block, good, 11, 21, 1; x0=10, y0=20) == (5.0, 0.0)
    @test isnan(first(FY4B._window_value(block, good, 10, 20, 3; x0=10, y0=20)))  # 4 of 9 inside
end

@testset "the shifted grid round-trips through its chunked NetCDF" begin
    mktempdir() do dir
        path = joinpath(dir, "fy4b_offset_grid.nc")
        offsets = (window=3, radius=1)
        ids = ["a", "b"]
        # More hours than one chunk, so a block boundary is crossed and the tail is flushed.
        total = FY4B.OFFSET_GRID_BLOCK + 5
        times = DateTime(2022, 6, 1, 9) .+ Hour.(0:(total - 1))
        hour_grid(h) = [h + 10i + 100a + 1000b for i in 1:2, a in 1:3, b in 1:3]   # [station, dx, dy]
        on_shifted, finish, _ = FY4B._offset_grid_writer(path, offsets, ids, times)
        for h in 1:total
            grid = Float64.(hour_grid(h))
            h == 7 && (grid[2, 3, 1] = NaN)                                  # a gap survives
            on_shifted(h, grid)
        end
        finish()
        @test isfile(path) && !any(startswith("fy4b_offset_grid.nc.tmp"), readdir(dir))
        NCDataset(path) do ds
            @test ds["offset_x"][:] == [-1, 0, 1] && ds["offset_y"][:] == [-1, 0, 1]
            @test ds["station_id"][:] == ids
            @test ds["time_label"][1] == "2022-06-01T09:00:00Z"
            grid = ds["fy4b"]
            @test size(grid) == (3, 3, 2, total)
            # [dx, dy, station, time] reads back the hour's [station, dx, dy] values.
            @test grid[3, 2, 1, total] == Float32(total + 10 + 300 + 2000)
            @test grid[1, 1, 2, FY4B.OFFSET_GRID_BLOCK] == Float32(FY4B.OFFSET_GRID_BLOCK + 20 + 100 + 1000)
            @test isnan(grid[3, 1, 2, 7]) && !isnan(grid[3, 1, 1, 7])
        end
    end
end

@testset "a broken copy sorting first does not shadow a good one" begin
    mktempdir() do dir
        name = "FY4B-_AGRI--_N_DISK_1050E_L2-_QPE-_MULT_NOM_20240608064500_20240608065959_4000M_V0001.NC"
        mkpath(joinpath(dir, "a")); mkpath(joinpath(dir, "b"))
        write(joinpath(dir, "a", name), "truncated download")             # sorts first, will not open
        NCDataset(joinpath(dir, "b", name), "c") do ds
            defDim(ds, "x", FY4B.GRID_SIZE); defDim(ds, "y", FY4B.GRID_SIZE)
            defVar(ds, "Precipitation", Float32, ("x", "y"); deflatelevel=1)
        end
        @test FY4B.segment_problem(joinpath(dir, "a", name)) isa String
        @test FY4B.segment_problem(joinpath(dir, "b", name)) === nothing
        groups = FY4B.find_nc_files(dir; years=2024:2024, months=6:6)
        @test only(groups[DateTime(2024, 6, 8, 6)]) == joinpath(dir, "b", name)
        # With no good copy the first is still kept, so the hour is reported unreadable, not lost.
        rm(joinpath(dir, "b", name))
        write(joinpath(dir, "b", name), "also broken")
        groups = FY4B.find_nc_files(dir; years=2024:2024, months=6:6)
        @test only(groups[DateTime(2024, 6, 8, 6)]) == joinpath(dir, "a", name)
        @test haskey(FY4B.find_unreadable_files(groups), joinpath(dir, "a", name))
    end
end
