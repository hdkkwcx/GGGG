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

@testset "an adaptive map is applied only as strongly as the training gauges support" begin
    rng = MersenneTwister(5)
    n_station, n_time = 10, 3000
    truth = [rand(rng) < 0.2 ? 4 * rand(rng)^2 : 0.0 for _ in 1:n_station, _ in 1:n_time]
    train = 1:7
    full = length(QM_STRENGTHS)
    # Twice the amounts and nothing else: the full map undoes it exactly, so it wins.
    scaled = 2 .* truth
    cells = fill(3, n_time)
    maps = fit_quantile_map(truth, scaled, train, cells; min_wet=100)
    choice = choose_strength(truth, scaled, maps, train, cells)
    @test choice[3].candidate == full
    @test choice[3].rmse_chosen == choice[3].rmse_full < choice[3].rmse_none
    # Small amounts unrelated to the gauges: restoring the gauge tail at random places only adds
    # error, so a smaller adjustment wins.
    unrelated = [rand(rng) < 0.2 ? 0.5 * rand(rng) : 0.0 for _ in 1:n_station, _ in 1:n_time]
    maps = fit_quantile_map(truth, unrelated, train, cells; min_wet=100)
    choice = choose_strength(truth, unrelated, maps, train, cells)
    @test choice[3].candidate < full
    @test choice[3].rmse_chosen < choice[3].rmse_full
    # Held-out gauges cannot move the choice.
    scrambled = copy(truth)
    scrambled[8:10, :] .= 100.0
    @test choose_strength(scrambled, unrelated, maps, train, cells)[3].candidate == choice[3].candidate
    # The choice can be finer than the map: one season map, two months. In the first the
    # satellite is the gauges' amounts doubled and the map helps; in the second it is unrelated
    # rain and the map hurts. A per-month choice keeps the first and holds back on the second.
    half = n_time ÷ 2
    mixed = hcat(scaled[:, 1:half], unrelated[:, (half + 1):end])
    months = [t <= half ? 1 : 2 for t in 1:n_time]
    season_maps = fit_quantile_map(truth, mixed, train, cells; min_wet=100)
    by_month = choose_strength(truth, mixed, season_maps, train, cells, months)
    @test sort(collect(keys(by_month))) == [1, 2]
    @test by_month[1].candidate == full
    @test by_month[2].candidate < full
    @test by_month[1].n_pairs + by_month[2].n_pairs == length(train) * n_time
    # Every candidate is monotone in the satellite value.
    m = maps[3]
    xs = collect(0.0:0.01:12.0)
    for c in 0:full
        @test issorted([SatelliteCalibration._adjusted_value(x, SatelliteCalibration._map_value(m, x),
            m.tau, c) for x in xs])
    end
end

@testset "a fold's anchor uses each regime's chosen shift, then the map" begin
    mktempdir() do dir
        s = synthetic_calibration(dir)
        grid = open_offset_grid(s.path, s.ids, s.times)
        labels, regime_of_time = product_regimes(s.meta, s.times)
        sums = shift_sums(grid, s.rain, regime_of_time, 2)
        cells = calibration_cells(s.times, regime_of_time)
        cache = Dict{Int,Matrix{Float64}}()
        anchor, shift_rows, map_rows, strength_rows = calibrate_fold(
            grid, sums, s.rain, 1:6, labels, regime_of_time, cells; cache)
        @test size(anchor) == size(s.rain)
        @test [(row.dx, row.dy) for row in shift_rows] == collect(s.true_shift)
        @test length(cache) == 2
        @test cor(vec(anchor), vec(s.rain)) > 0.9
        @test Set(row.regime for row in map_rows) == Set(labels)
        # Without `adaptive`, every cell gets the full map: the anchor is the plain map's.
        shifted = fill(NaN, size(s.rain))
        for (g, row) in enumerate(shift_rows)
            k = findfirst(==((row.dx, row.dy)), grid.shifts)
            cols = findall(==(g), regime_of_time)
            shifted[:, cols] = read_shift(grid, k)[:, cols]
        end
        plain = apply_quantile_map(fit_quantile_map(s.rain, shifted, 1:6, cells), shifted, cells)
        @test isequal(anchor, plain)
        @test !hasproperty(first(map_rows), :strength)
        @test all(row.strength == 1.0 && row.threshold_applied for row in strength_rows)
        @test all(row.train_rmse_chosen == row.train_rmse_full for row in strength_rows)
        # Adaptive, with the strength chosen per month cell as the benchmark does.
        month_labels, months = strength_cells(s.times, regime_of_time)
        @test month_labels == ["2022-06", "2022-06"]           # one month, split by the regimes
        adapted, _, adapted_maps, adapted_rows = calibrate_fold(
            grid, sums, s.rain, 1:6, labels, regime_of_time, cells; cache, adaptive=true,
            choice_of_time=months, choice_labels=month_labels)
        @test size(adapted) == size(anchor)
        @test adapted_maps == map_rows
        @test [(row.regime, row.cell) for row in adapted_rows] ==
            [(labels[1], "2022-06"), (labels[2], "2022-06")]
        @test all(row.strength in QM_STRENGTHS for row in adapted_rows)
        @test all(row.train_rmse_chosen <= min(row.train_rmse_none, row.train_rmse_full)
                  for row in adapted_rows)
    end
end

@testset "a crossed anchor never uses the gauges of its own weeks" begin
    mktempdir() do dir
        s = synthetic_calibration(dir)
        grid = open_offset_grid(s.path, s.ids, s.times)
        labels, regime_of_time = product_regimes(s.meta, s.times)
        cells = calibration_cells(s.times, regime_of_time)
        month_labels, months = strength_cells(s.times, regime_of_time)
        # Each regime's first half is time fold 1 and its second half fold 2, so that after the
        # guard days both folds still train on both regimes.
        fold_of_time = [(h - 1) % 200 < 100 ? 1 : 2 for h in eachindex(s.times)]
        for fold in 1:2
            kept = Main.TimeBlockFolds.training_hours(fold_of_time, s.times, fold)
            @test all(count(kept .& (regime_of_time .== g)) >= 24 for g in 1:2)
        end
        train = 1:6
        crossed(y) = calibrate_crossed(grid,
            crossed_shift_sums(grid, y, regime_of_time, 2, fold_of_time, s.times), y, train,
            labels, regime_of_time, cells, fold_of_time, s.times;
            adaptive=true, choice_of_time=months, choice_labels=month_labels)
        anchor, shift_rows, map_rows, strength_rows = crossed(s.rain)
        @test size(anchor) == size(s.rain) && all(isfinite, anchor)
        @test Set(row.time_fold for row in map_rows) == Set([1, 2])
        # Each time fold still finds every regime's planted shift, on its own training hours.
        @test length(shift_rows) == 2 * length(labels)
        @test all((row.dx, row.dy) == s.true_shift[findfirst(==(row.regime), labels)]
                  for row in shift_rows)
        for fold in 1:2
            kept = Main.TimeBlockFolds.training_hours(fold_of_time, s.times, fold)
            hours = findall(==(fold), fold_of_time)
            @test !any(kept[hours])
            # The monthly strength was chosen on the training hours alone.
            @test sum(row.n_pairs for row in strength_rows if row.time_fold == fold) ==
                length(train) * count(kept)
            # Rain anywhere in a fold's own hours or its guard days cannot move its anchor. It
            # does move the other fold's, whose training hours they are.
            noisy = copy(s.rain)
            noisy[:, .!kept] .= 50.0
            moved = first(crossed(noisy))
            @test isequal(moved[:, hours], anchor[:, hours])
            others = findall(!=(fold), fold_of_time)
            @test !isequal(moved[:, others], anchor[:, others])
        end
    end
end
