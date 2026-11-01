using Test, CSV, DataFrames, Dates, Random, Statistics

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("FY4BPreprocessing", "SatelliteTwoStage")
using Main.SatelliteCalibration
using Main.SatelliteTwoStage

# An offset grid of radius 1 (shifts -1:1 in x and y) from `value(station, hour, dx, dy)`.
function two_stage_grid(dir, ids, times, value)
    path = joinpath(dir, "fy4b_offset_grid.nc")
    on_shifted, finish, _ = Main.FY4BPreprocessing._offset_grid_writer(
        path, (window=3, radius=1), ids, times)
    for h in eachindex(times)
        on_shifted(h, [value(i, h, dx, dy) for i in eachindex(ids), dx in -1:1, dy in -1:1])
    end
    finish()
    return open_offset_grid(path, ids, times)
end

@testset "context classes read intensity, rain nearby and the adjacent hours" begin
    @test fy4b_intensity_class.([0.0, 0.099, 0.1, 0.5, 1.0, 15.9, 16.0, 100.0]) == [1, 1, 2, 3, 4, 7, 8, 8]
    @test context_label(1) == (; intensity="dry", nearby="none", adjacent="neither")
    @test context_label(N_CONTEXT) == (; intensity=">=16", nearby="all", adjacent="both")
    @test context_label(context_class(4, 2, 1)) == (; intensity="1-2", nearby="some", adjacent="one")
    @test N_CONTEXT == 72 && length(unique(context_class(i, j, m)
        for i in 1:8, j in 1:3, m in 0:2)) == 72
    # Hours 0, 1, 2, then a gap to 5 and 6: the gap's neighbours are not adjacent.
    times = DateTime(2022, 6, 1) .+ Hour.([0, 1, 2, 5, 6])
    shifted = [0.3 2.5 0.0 1.0 0.05
               0.0 NaN 0.0 0.0 0.0]
    nearby = UInt8[1 2 3 1 2
                   1 1 1 1 1]
    classes = context_classes(shifted, nearby, times)
    @test classes[1, :] == [context_class(2, 1, 1), context_class(5, 2, 1),
        context_class(1, 3, 1), context_class(4, 1, 0), context_class(1, 2, 1)]
    # A missing value has no class and counts as dry for its neighbours.
    @test classes[2, :] == [context_class(1, 1, 0), 0, context_class(1, 1, 0),
        context_class(1, 1, 0), context_class(1, 1, 0)]
    @test_throws ArgumentError context_classes(shifted, nearby, reverse(times))
end

@testset "rain nearby counts the wet shifts around the chosen one" begin
    mktempdir() do dir
        ids = ["a", "b", "c"]
        times = DateTime(2022, 6, 1) .+ Hour.(0:3)
        # a: wet everywhere except one gap in the last hour; b: wet only at (1, 1); c: dry.
        value(i, h, dx, dy) = i == 1 ? ((h, dx, dy) == (4, -1, -1) ? NaN : 1.0) :
            i == 2 ? ((dx, dy) == (1, 1) ? 1.0 : 0.0) : 0.0
        grid = two_stage_grid(dir, ids, times, value)
        centre = findfirst(==((0, 0)), grid.shifts)
        east = findfirst(==((1, 0)), grid.shifts)
        count, n_neighbour = neighbour_wet_count(grid, centre)
        @test n_neighbour == 8
        @test count[1, :] == [8, 8, 8, 7] && all(==(1), count[2, :]) && all(==(0), count[3, :])
        # Shifts past the grid's edge are not neighbours: (1, 0) has five inside it.
        @test neighbour_wet_count(grid, east)[2] == 5
        cache = Dict{Int,Tuple{Matrix{UInt8},Int}}()
        nearby = fold_nearby(grid, [centre, east, 0], [1, 1, 2, 2]; cache)
        @test nearby[:, 1:2] == UInt8[3 3; 2 2; 1 1]
        @test nearby[:, 3:4] == UInt8[3 3; 2 2; 1 1]      # all five of (1, 0)'s; (1, 1) is one
        @test length(cache) == 2
        @test fold_nearby(grid, [centre, 0], [1, 2, 2, 2])[:, 2:4] == zeros(UInt8, 3, 3)
    end
end

@testset "the tables are counts on the training stations, pooled upward when sparse" begin
    n_station, n_time = 4, 100
    train = 1:3
    # JJA and SON of one regime, then JJA of a second.
    cells = [t <= 60 ? 3 : t <= 90 ? 4 : 7 for t in 1:n_time]
    q_a = context_class(5, 3, 2)               # 2-4 mm/h, rain all around and in both hours
    q_b = context_class(1, 1, 0)               # dry and isolated
    q_rare = context_class(5, 1, 0)            # once, in the first cell
    q_absent = context_class(8, 1, 0)          # never
    classes = [iseven(t) ? q_a : q_b for _ in 1:n_station, t in 1:n_time]
    classes[1, 1] = q_rare
    classes[2, 3] = 0                          # no FY4B value
    y = [classes[s, t] in (q_a, q_rare) ? (s <= 2 ? 3.0 : 0.0) :
         (s == 1 && t % 6 == 1 ? 0.2 : 0.0) for s in 1:n_station, t in 1:n_time]
    y[3, 5] = NaN                              # no gauge value
    tables = fit_context_tables(y, classes, cells, train; min_pairs=50)
    @test tables.cells == [3, 4, 7]
    own(cell, q) = count(classes[s, t] == q && cells[t] == cell && !isnan(y[s, t])
        for s in train, t in 1:n_time)
    own_wet(cell, q) = count(classes[s, t] == q && cells[t] == cell && y[s, t] >= 0.1
        for s in train, t in 1:n_time)
    for cell in (3, 4, 7), q in (q_a, q_b)
        @test tables.n_own[cell, q] == own(cell, q)
        @test tables.n_wet_own[cell, q] == own_wet(cell, q)
    end
    pool(cell, q) = POOL_LEVELS[tables.pool[cell, q]]
    @test pool(3, q_a) == "cell"                                 # 90 pairs of its own
    @test pool(4, q_a) == "regime"                               # 45, but 135 in its regime
    @test pool(7, q_a) == "all"                                  # 15, and its regime has no more
    @test pool(3, q_rare) == "intensity" && pool(3, q_absent) == "pooled"
    @test tables.p[3, q_a] == 60 / 90 && tables.amount[3, q_a] == 3.0
    @test tables.n[4, q_a] == own(3, q_a) + own(4, q_a)
    # P x amount is the class's mean gauge rain on the training stations (dry hours are 0 here).
    expected = expected_anchor(tables, classes, cells)
    @test expected[4, 2] ≈ mean(y[s, t] for s in train, t in 1:60 if classes[s, t] == q_a)
    @test expected[1, 1] == tables.p[3, q_rare] * tables.amount[3, q_rare]
    @test isnan(expected[2, 3])
    # The occurrence rule: in the first cell the classes with the highest P cover the gauges' 70
    # wet pairs once q_a is in, so q_rare and q_a are wet and q_b is not.
    target = sum(tables.n_wet_own[3, :])
    @test target == 70
    @test tables.p_star[3] == tables.p[3, q_a] > tables.p[3, q_b]
    wet_pairs = tables.n_own[3, q_rare] + tables.n_own[3, q_a]
    @test wet_pairs >= target > wet_pairs - tables.n_own[3, q_a]
    calibrated = [s == 1 ? 5.0 : s == 2 ? 0.0 : 0.4 for s in 1:n_station, t in 1:n_time]
    calibrated[4, 4] = NaN
    occurrence = occurrence_anchor(tables, classes, cells, calibrated)
    @test occurrence[1, 2] == 5.0                                # wet: keeps the calibrated amount
    @test occurrence[2, 2] == tables.amount[3, q_a]              # wet, calibrated 0: the class mean
    @test occurrence[3, 1] == 0.0 && occurrence[1, 3] == 0.0     # q_b is not called wet
    @test isnan(occurrence[2, 3]) && isnan(occurrence[4, 4])
    # Held-out gauges cannot move the tables or either anchor.
    scrambled = copy(y)
    scrambled[4, :] .= 100.0
    again = fit_context_tables(scrambled, classes, cells, train; min_pairs=50)
    @test isequal(expected_anchor(again, classes, cells), expected)
    @test isequal(occurrence_anchor(again, classes, cells, calibrated), occurrence)
    rows = context_table_rows(tables, ["133.0E_V1.0", "105.0E_V1.0.0"])
    @test length(rows) == 3 * N_CONTEXT
    @test first(rows).regime == "133.0E_V1.0" && first(rows).season == "JJA"
    @test rows[end].regime == "105.0E_V1.0.0"
    @test count(row -> row.called_wet && row.season == "JJA" && row.regime == "133.0E_V1.0",
        rows) >= 2
    @test_throws ErrorException fit_context_tables(y, zeros(Int, size(y)), cells, train)
    # A class that overshoots the wet count by more than it is short stays dry: 10 likely-wet
    # hours cover 10 of the 25 wet pairs, and adding 990 rarely-wet ones would cover 1000.
    q_hi, q_lo = context_class(4, 3, 2), context_class(1, 1, 0)
    sparse_classes = reshape([t <= 10 ? q_hi : q_lo for t in 1:1000], 1, :)
    sparse_y = reshape([t <= 10 ? (isodd(t) ? 2.0 : 0.0) : (t % 50 == 0 ? 1.0 : 0.0)
        for t in 1:1000], 1, :)
    sparse = fit_context_tables(sparse_y, sparse_classes, fill(3, 1000), [1]; min_pairs=1)
    @test sum(sparse.n_wet_own[3, :]) == 25
    @test sparse.p_star[3] == sparse.p[3, q_hi] == 0.5
    @test occurrence_anchor(sparse, sparse_classes, fill(3, 1000), fill(0.0, 1, 1000))[11:end] ==
        zeros(990)
end

@testset "a fold's two-stage anchors come from the fold's shift and training gauges" begin
    mktempdir() do dir
        rng = MersenneTwister(11)
        n_station, n_time = 8, 400
        ids = ["s$i" for i in 1:n_station]
        times = DateTime(2022, 6, 1, 9) .+ Hour.(0:(n_time - 1))
        rain = [rand(rng) < 0.3 ? 5 * rand(rng) : 0.0 for _ in 1:n_station, _ in 1:n_time]
        true_shift = ((1, 0), (-1, 1))
        regime_of(h) = h <= n_time ÷ 2 ? 1 : 2
        noise = [rand(rng) < 0.3 ? 5 * rand(rng) : 0.0 for _ in 1:n_station, _ in 1:n_time, _ in 1:9]
        value(i, h, dx, dy) = (dx, dy) == true_shift[regime_of(h)] ? rain[i, h] + 0.1 :
            noise[i, h, (dy + 1) * 3 + dx + 2]
        grid = two_stage_grid(dir, ids, times, value)
        meta = joinpath(dir, "fy4b_hour_meta.csv")
        CSV.write(meta, DataFrame(
            time=Dates.format.(times, dateformat"yyyy-mm-ddTHH:MM:SS") .* "Z",
            sat_lon=fill(133.0, n_time),
            software=[regime_of(h) == 1 ? "V1.0" : "V1.0.0" for h in 1:n_time]))
        labels, regime_of_time = product_regimes(meta, times)
        sums = shift_sums(grid, rain, regime_of_time, 2)
        cells = calibration_cells(times, regime_of_time)
        train = 1:6
        cache = Dict{Int,Matrix{Float64}}()
        # `shift_fold` is what `calibrate_fold` shifts with.
        shifted, chosen, shift_rows = shift_fold(grid, sums, train, labels, regime_of_time; cache)
        @test [grid.shifts[k] for k in chosen] == collect(true_shift)
        @test [(row.dx, row.dy) for row in shift_rows] == collect(true_shift)
        for g in 1:2
            cols = findall(==(g), regime_of_time)
            @test isequal(shifted[:, cols], read_shift(grid, chosen[g])[:, cols])
        end
        calibrated = first(calibrate_fold(grid, sums, rain, train, labels, regime_of_time, cells;
            cache))
        nearby_cache = Dict{Int,Tuple{Matrix{UInt8},Int}}()
        expected, occurrence, rows = two_stage_fold(grid, sums, rain, train, labels,
            regime_of_time, cells, times, calibrated; cache, nearby_cache, min_pairs=20)
        @test size(expected) == size(occurrence) == size(rain)
        @test all(isfinite, expected) && all(isfinite, occurrence)
        @test length(nearby_cache) == 2
        @test length(rows) == length(unique(cells)) * N_CONTEXT
        @test cor(vec(expected), vec(rain)) > 0.9
        # Held-out gauges cannot move either anchor.
        scrambled = copy(rain)
        scrambled[7:8, :] .= reverse(scrambled[7:8, :]; dims=2)
        again = two_stage_fold(grid, sums, scrambled, train, labels, regime_of_time, cells, times,
            calibrated; cache, nearby_cache, min_pairs=20)
        @test isequal(again[1], expected) && isequal(again[2], occurrence)
    end
end
