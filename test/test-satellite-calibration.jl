using Test, CSV, DataFrames, Dates, Random, Statistics

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("FY4BPreprocessing", "SatelliteCalibration")
using Main.SatelliteCalibration
const FY4B_SOURCE = Main.FY4BPreprocessing

# A synthetic grid where the right shift is known: in the first regime the gauge field appears at
# (+1, 0), in the second at (-1, +1); every other shift holds unrelated rain.
function synthetic_calibration(dir; n_station=8, n_time=400)
    rng = MersenneTwister(7)
    times = DateTime(2022, 6, 1, 9) .+ Hour.(0:(n_time - 1))
    rain = [rand(rng) < 0.3 ? 5 * rand(rng) : 0.0 for _ in 1:n_station, _ in 1:n_time]
    regime_of(h) = h <= n_time ÷ 2 ? 1 : 2
    true_shift = ((1, 0), (-1, 1))
    ids = ["s$i" for i in 1:n_station]
    path = joinpath(dir, "fy4b_offset_grid.nc")
    on_shifted, finish, _ = FY4B_SOURCE._offset_grid_writer(path, (window=3, radius=1), ids, times)
    for h in 1:n_time
        grid = [(dx, dy) == true_shift[regime_of(h)] ? rain[i, h] + 0.1 * rand(rng) :
                (rand(rng) < 0.3 ? 5 * rand(rng) : 0.0)
                for i in 1:n_station, dx in -1:1, dy in -1:1]
        on_shifted(h, grid)
    end
    finish()
    meta = joinpath(dir, "fy4b_hour_meta.csv")
    CSV.write(meta, DataFrame(
        time=Dates.format.(times, dateformat"yyyy-mm-ddTHH:MM:SS") .* "Z",
        sat_lon=fill(133.0, n_time),
        software=[regime_of(h) == 1 ? "V1.0" : "V1.0.0" for h in 1:n_time],
    ))
    return (; times, rain, ids, path, meta, true_shift)
end

@testset "regimes come from the hour metadata, not from dates" begin
    mktempdir() do dir
        s = synthetic_calibration(dir)
        labels, regime_of_time = product_regimes(s.meta, s.times)
        @test labels == ["133.0E_V1.0", "133.0E_V1.0.0"]
        @test regime_of_time == [h <= 200 ? 1 : 2 for h in 1:400]
        @test_throws ErrorException product_regimes(s.meta, [s.times; s.times[end] + Hour(1)])
        cells = calibration_cells(s.times, regime_of_time)
        @test all(==(3), cells[1:200]) && all(==(7), cells[201:end])      # JJA in both regimes
    end
end

@testset "the shift is chosen per regime on training stations only" begin
    mktempdir() do dir
        s = synthetic_calibration(dir)
        grid = open_offset_grid(s.path, s.ids, s.times)
        @test length(grid.shifts) == 9 && grid.shifts[5] == (0, 0)
        _, regime_of_time = product_regimes(s.meta, s.times)
        sums = shift_sums(grid, s.rain, regime_of_time, 2)
        train = 1:6
        for g in 1:2
            k, r, r0 = choose_shift(grid, sums, train, g)
            @test grid.shifts[k] == s.true_shift[g]
            @test r > 0.9 > r0
        end
        # Held-out gauges cannot move the choice: scramble them and choose again.
        scrambled = copy(s.rain)
        scrambled[7:8, :] .= reverse(scrambled[7:8, :]; dims=2)
        rescored = shift_sums(grid, scrambled, regime_of_time, 2)
        @test all(choose_shift(grid, rescored, train, g)[1] == choose_shift(grid, sums, train, g)[1]
                  for g in 1:2)
        @test read_shift(grid, 5) == read_shift(open_offset_grid(s.path, s.ids, s.times), 5)
    end
end

@testset "the quantile map matches the training gauges' wet frequency and amounts" begin
    rng = MersenneTwister(3)
    n_station, n_time = 10, 3000
    truth = [rand(rng) < 0.2 ? 4 * rand(rng)^2 : 0.0 for _ in 1:n_station, _ in 1:n_time]
    # Twice the amounts, and false rain in a quarter of the dry hours.
    satellite = [v > 0 ? 2v : (rand(rng) < 0.25 ? 0.3 * rand(rng) : 0.0) for v in truth]
    cells = fill(3, n_time)
    train = 1:7
    maps = fit_quantile_map(truth, satellite, train, cells; min_wet=100)
    mapped = apply_quantile_map(maps, satellite, cells)
    @test maps[3].pool == "cell"
    @test isapprox(mean(mapped[train, :] .>= 0.1), mean(truth[train, :] .>= 0.1); atol=0.01)
    @test isapprox(mean(mapped[train, :]), mean(truth[train, :]); rtol=0.05)
    @test mean(satellite[train, :]) > 1.5 * mean(truth[train, :])            # it did need fixing
    # Monotone: a larger satellite value never maps lower.
    xs = collect(0.0:0.01:12.0)
    ys = [SatelliteCalibration._map_value(maps[3], x) for x in xs]
    @test issorted(ys)
    @test ys[end] > maps[3].qy[end]                                           # the cap stays a lower bound
    # Held-out gauges cannot move the map.
    scrambled = copy(truth)
    scrambled[8:10, :] .= 100.0
    @test fit_quantile_map(scrambled, satellite, train, cells; min_wet=100)[3].qy == maps[3].qy
    # A cell too sparse for its own map takes its regime's.
    sparse = fit_quantile_map(truth, satellite, train, [t <= 50 ? 4 : 3 for t in 1:n_time];
        min_wet=100)
    @test sparse[4].pool == "regime" && sparse[3].pool == "cell"
end

@testset "a fold's anchor uses each regime's chosen shift, then the map" begin
    mktempdir() do dir
        s = synthetic_calibration(dir)
        grid = open_offset_grid(s.path, s.ids, s.times)
        labels, regime_of_time = product_regimes(s.meta, s.times)
        sums = shift_sums(grid, s.rain, regime_of_time, 2)
        cells = calibration_cells(s.times, regime_of_time)
        cache = Dict{Int,Matrix{Float64}}()
        anchor, shift_rows, map_rows = calibrate_fold(
            grid, sums, s.rain, 1:6, labels, regime_of_time, cells; cache)
        @test size(anchor) == size(s.rain)
        @test [(row.dx, row.dy) for row in shift_rows] == collect(s.true_shift)
        @test length(cache) == 2
        @test cor(vec(anchor), vec(s.rain)) > 0.9
        @test Set(row.regime for row in map_rows) == Set(labels)
    end
end
