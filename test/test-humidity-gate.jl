using Test, Dates, Random

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("HumidityGate", "TimeBlockFolds")
const HG = Main.HumidityGate
const HG_TBF = Main.TimeBlockFolds

"""
A network where rain reaches the gauges only in moist air (RH >= 75), and the anchor also rains
lightly in dry air (RH < 60), plus a few heavy dry-air values the gate must leave alone.
"""
function hg_fixture(; n_station=20, n_day=60, seed=3)
    rng = MersenneTwister(seed)
    times = [DateTime(2023, 6, 1, 1) + Hour(h) for h in 0:(24 * n_day - 1)]
    rh = 40.0 .+ 60.0 .* rand(rng, n_station, length(times))
    y = [r >= 75 && rand(rng) < 0.5 ? 1.0 + 3rand(rng) : 0.0 for r in rh]
    x = [yv > 0 ? yv : (r < 60 && rand(rng) < 0.4 ? 0.5 + rand(rng) : 0.0) for (yv, r) in zip(y, rh)]
    heavy = findall(r -> r < 50, rh)[1:5]
    x[heavy] .= 12.0
    return (; times, rh, y, x, heavy)
end

@testset "HumidityGate" begin
    f = hg_fixture()
    stations, hours = collect(1:15), collect(axes(f.x, 2))

    @testset "the cut removes dry-air false rain and spares heavy values" begin
        choice = HG.choose_gate_cut(f.x, f.y, f.rh, stations, hours)
        @test 60 <= choice.cut <= 75
        @test choice.csi_after > choice.csi_before
        gated = HG.apply_gate(f.x, f.rh, choice.cut)
        @test all(gated[f.heavy] .== 12.0)
        dry = (f.rh .< 60) .& (f.x .< 8)
        @test all(gated[dry] .== 0.0)
        @test all(gated[f.rh .>= choice.cut] .== f.x[f.rh .>= choice.cut])
    end

    @testset "no gate when no cut helps" begin
        honest = copy(f.y)                      # wet only where the gauge is: no false rain
        choice = HG.choose_gate_cut(honest, f.y, f.rh, stations, hours)
        @test choice.cut == -Inf && choice.csi_after == choice.csi_before
        @test HG.apply_gate(honest, f.rh, choice.cut) == honest
    end

    @testset "NaN never gates and is never counted" begin
        x, rh = copy(f.x), copy(f.rh)
        x[1, 1], rh[2, 2] = NaN, NaN
        gated = HG.apply_gate(x, rh, 90.0)
        @test isnan(gated[1, 1]) && gated[2, 2] == x[2, 2]
        @test HG.choose_gate_cut(x, f.y, rh, [1, 2], hours).n_cells == 2 * length(hours) - 2
    end

    @testset "the cut sees training stations and training weeks only" begin
        folds = HG_TBF.time_block_folds(f.times, fill("all", length(f.times)); seed=7)
        gated, rows = HG.time_crossed_gate(f.x, f.y, f.rh, stations, folds.fold_of_time, f.times)
        @test length(rows) == 5 && all(row.n_train_hours < length(f.times) for row in rows)
        # Held-out stations: changing their gauges and anchor changes no cut.
        y2, x2 = copy(f.y), copy(f.x)
        y2[16:end, :] .= 50.0
        x2[16:end, :] .= 0.0
        _, rows2 = HG.time_crossed_gate(x2, y2, f.rh, stations, folds.fold_of_time, f.times)
        @test [row.cut for row in rows] == [row.cut for row in rows2]
        # Held-out weeks: wrecking fold 1's hours leaves fold 1's cut alone.
        y3 = copy(f.y)
        y3[:, folds.fold_of_time .== 1] .= 50.0
        _, rows3 = HG.time_crossed_gate(f.x, y3, f.rh, stations, folds.fold_of_time, f.times)
        @test rows3[1].cut == rows[1].cut
        # Each fold's hours carry that fold's own cut, at every station.
        for row in rows
            columns = folds.fold_of_time .== row.time_fold
            @test gated[:, columns] == HG.apply_gate(f.x[:, columns], f.rh[:, columns], row.cut)
        end
        @test HG.time_crossed_gate(f.x, f.y, f.rh, stations, folds.fold_of_time, f.times)[1] == gated
    end

    @testset "station_humidity is ERA5LandProcessing's RH" begin
        t2m, d2m = [20.0 25.0], [20.0 15.0]
        rh = HG.station_humidity(t2m, d2m)
        @test rh[1] ≈ 100.0 && 50 < rh[2] < 60
    end
end
