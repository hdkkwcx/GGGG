using Test, Dates, Random, Statistics

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("TimeBlockFolds")
const TBF = Main.TimeBlockFolds

"""Hours of `n_days` days, labelled by their ends, with no hour missing."""
tbf_hours(n_days) = [DateTime(2022, 6, 1, 1) + Hour(h) for h in 0:(24 * n_days - 1)]

"""
Accuracy on held-out cells of a table of P(wet) per class, fitted on each split's training cells.
A class the training cells never saw takes their overall rate.
"""
function tbf_table_accuracy(class, wet, splits)
    correct = total = 0
    for split in splits
        counts, wets = Dict{Int,Int}(), Dict{Int,Int}()
        for t in split.train_hours, i in split.train_stations
            counts[class[i, t]] = get(counts, class[i, t], 0) + 1
            wets[class[i, t]] = get(wets, class[i, t], 0) + wet[i, t]
        end
        overall = sum(values(wets)) / sum(values(counts))
        for t in split.test_hours, i in split.test_stations
            c = class[i, t]
            p = haskey(counts, c) ? wets[c] / counts[c] : overall
            correct += (p >= 0.5) == wet[i, t]
            total += 1
        end
    end
    return correct / total
end

@testset "TimeBlockFolds" begin
    times = tbf_hours(400)
    strata = [month(TBF.label_day(t)) <= 8 ? "summer" : "rest" for t in times]
    folds = TBF.time_block_folds(times, strata; seed=7)

    @testset "blocks are whole, and a seed fixes the deal" begin
        for block in unique(folds.block_of_time)
            hours = findall(==(block), folds.block_of_time)
            @test length(unique(TBF.label_day.(times[hours]))) <= 7
            @test length(unique(folds.fold_of_time[hours])) == 1
        end
        @test folds.fold_of_time == TBF.time_block_folds(times, strata; seed=7).fold_of_time
        @test folds.fold_of_time != TBF.time_block_folds(times, strata; seed=8).fold_of_time
        @test sort(unique(folds.fold_of_time)) == 1:5
    end

    @testset "every stratum spreads evenly over the folds" begin
        for stratum in ("summer", "rest")
            chosen = strata .== stratum
            sizes = [length(unique(folds.block_of_time[chosen .& (folds.fold_of_time .== f)]))
                for f in 1:5]
            @test maximum(sizes) - minimum(sizes) <= 2
            @test minimum(sizes) >= 1
        end
        sizes = [length(unique(folds.block_of_time[folds.fold_of_time .== f])) for f in 1:5]
        @test maximum(sizes) - minimum(sizes) <= 1
    end

    @testset "a gap in the record does not leave an empty block" begin
        gappy = vcat(times[1:240], times[(24 * 100 + 1):end])
        g = TBF.time_block_folds(gappy, fill("all", length(gappy)); seed=7)
        @test sort(unique(g.block_of_time)) == 1:maximum(g.block_of_time)
        @test length(g.fold_of_block) == maximum(g.block_of_time)
    end

    @testset "training stops a guard short of held-out days" begin
        for fold in 1:5, guard_days in (0, 1, 3)
            train = TBF.training_hours(folds.fold_of_time, times, fold; guard_days)
            held = TBF.label_day.(times[folds.fold_of_time .== fold])
            for day in TBF.label_day.(times[train])
                @test minimum(abs(Dates.value(day - other)) for other in unique(held)) > guard_days
            end
            @test !any(train .& (folds.fold_of_time .== fold))
        end
        # A training hour is lost only to the guard: with none, it is every other fold's hour.
        @test TBF.training_hours(folds.fold_of_time, times, 2; guard_days=0) == (folds.fold_of_time .!= 2)
    end

    fold_of_station = repeat(1:5, 8)
    splits = TBF.crossed_splits(fold_of_station, folds.fold_of_time, times)

    @testset "every gauge-hour is scored exactly once, by a model that saw neither" begin
        @test length(splits) == 25
        scored = zeros(Int, length(fold_of_station), length(times))
        for split in splits
            scored[split.test_stations, split.test_hours] .+= 1
            @test isempty(intersect(split.train_stations, split.test_stations))
            @test isempty(intersect(split.train_hours, split.test_hours))
            @test all(fold_of_station[split.test_stations] .== split.gauge_fold)
            @test all(folds.fold_of_time[split.test_hours] .== split.time_fold)
        end
        @test all(==(1), scored)
    end

    # What the design is for. Every gauge shares an episode's fate, and episodes last 2 days.
    n_station, n_hour = length(fold_of_station), length(times)
    episode_of_hour = [(Dates.value(TBF.label_day(t) - TBF.label_day(first(times))) ÷ 2) + 1 for t in times]
    fate = rand(MersenneTwister(11), Bool, maximum(episode_of_hour))
    gauges_only = [(; gauge_fold=s, time_fold=0, train_stations=findall(!=(s), fold_of_station),
        test_stations=findall(==(s), fold_of_station), train_hours=collect(1:n_hour),
        test_hours=collect(1:n_hour)) for s in 1:5]

    @testset "a table that memorises episodes passes gauges-only and fails crossed" begin
        class = [episode_of_hour[t] for _ in 1:n_station, t in 1:n_hour]
        wet = [fate[episode_of_hour[t]] for _ in 1:n_station, t in 1:n_hour]
        for split in splits
            @test isempty(intersect(episode_of_hour[split.train_hours], episode_of_hour[split.test_hours]))
        end
        @test tbf_table_accuracy(class, wet, gauges_only) == 1.0
        @test tbf_table_accuracy(class, wet, splits) < 0.65
    end

    @testset "block_times groups whole blocks of dates, and one day changes nothing" begin
        gappy = vcat(times[1:240], times[(24 * 100 + 1):end])
        @test Date.(TBF.block_times(gappy; block_days=1)) == Date.(gappy)
        weekly = TBF.block_times(gappy; block_days=7)
        for (t, start) in zip(gappy, weekly)
            @test Date(start) <= Date(t) < Date(start) + Day(7)
            @test Dates.value(Date(start) - Date(minimum(gappy))) % 7 == 0
        end
        @test issorted(weekly)
        @test_throws ArgumentError TBF.block_times(gappy; block_days=0)
        @test TBF.block_days_option(["run_dir"]) == (nothing, ["run_dir"])
        @test TBF.block_days_option(["--block-days", "1"]) == (1, String[])
        @test TBF.block_days_option(["--block-days", "7", "run_dir"]) == (7, ["run_dir"])
        @test TBF.block_days_option(["run_dir", "--block-days", "14"]) == (14, ["run_dir"])
        @test_throws ArgumentError TBF.block_days_option(["--block-days"])
        @test_throws ArgumentError TBF.block_days_option(["--block-days", "0"])
        @test TBF.from_option(["run_dir"]) == (nothing, ["run_dir"])
        @test TBF.from_option(["--from", "2022-06", "x"]) == (Date(2022, 6, 1), ["x"])
        @test TBF.from_suffix(Date(2022, 6, 1)) == "_from202206" && TBF.from_suffix(nothing) == ""
        @test_throws ArgumentError TBF.from_option(["--from", "2022-6"])
        @test_throws ArgumentError TBF.from_option(["--from"])
    end

    @testset "a signal that generalises survives both" begin
        rng = MersenneTwister(12)
        class = rand(rng, 1:2, n_station, n_hour)
        wet = [rand(rng) < (c == 1 ? 0.9 : 0.1) for c in class]
        gauges = tbf_table_accuracy(class, wet, gauges_only)
        crossed = tbf_table_accuracy(class, wet, splits)
        @test isapprox(gauges, 0.9; atol=0.01)
        @test isapprox(crossed, gauges; atol=0.01)
    end
end
