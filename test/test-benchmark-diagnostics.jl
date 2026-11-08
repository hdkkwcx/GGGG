using Test, CSV, DataFrames, Dates

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("BenchmarkDiagnostics")
const BD = Main.BenchmarkDiagnostics

"""
A fixture whose NaN pattern is known by construction.

4 stations x 6 hours, one fold per pair of stations. Every cell is evaluable (gauge and
satellite finite) except one, which is knocked out on the satellite side so the evaluable
denominator differs from the raw cell count. `mgwr` then drops:

- hour 2 entirely for fold 1 (both its stations) — a whole-hour dropout,
- hour 4 for a single station of fold 1 — a scattered dropout,
- nothing at all for fold 2.
"""
function dropout_fixture()
    ids = ["s1", "s2", "s3", "s4"]
    fold_map = Dict("s1" => 1, "s2" => 1, "s3" => 2, "s4" => 2)
    y_obs = fill(1.0, 4, 6)
    y_sat = fill(0.5, 4, 6)
    # Station 1, hour 6 is not evaluable: the satellite is missing there.
    y_sat[1, 6] = NaN
    mgwr = fill(0.7, 4, 6)
    mgwr[1:2, 2] .= NaN      # whole hour for fold 1
    mgwr[1, 4] = NaN         # scattered
    clean = fill(0.7, 4, 6)
    predictions = Dict{String,Matrix{Float64}}("mgwr" => mgwr, "clean" => clean)
    return (; ids, fold_map, y_obs, y_sat, predictions)
end

@testset "dropout_table splits whole-hour from scattered failures" begin
    fixture = dropout_fixture()
    table = BD.dropout_table(
        fixture.y_obs, fixture.y_sat, fixture.predictions, fixture.ids, fixture.fold_map;
        scheme="balanced_spatial", product="GSMaP", methods=["mgwr", "clean"],
    )
    @test nrow(table) == 4

    fold1 = only(filter(row -> row.method == "mgwr" && row.fold == 1, table))
    # Fold 1 has 2 stations x 6 hours = 12 cells, minus the one non-evaluable satellite cell.
    @test fold1.n_evaluable == 11
    @test fold1.n_dropped == 3
    @test fold1.cells_in_all_dropped_hours == 2
    @test fold1.cells_scattered == 1
    @test fold1.hours_evaluable == 6
    @test fold1.hours_all_dropped == 1
    @test fold1.own_coverage ≈ 1 - 3 / 11

    fold2 = only(filter(row -> row.method == "mgwr" && row.fold == 2, table))
    @test fold2.n_dropped == 0
    @test fold2.own_coverage == 1.0

    for row in eachrow(filter(:method => ==("clean"), table))
        @test row.n_dropped == 0
        @test row.hours_all_dropped == 0
        @test row.own_coverage == 1.0
    end
end

@testset "dropout_table counts an hour as fully dropped only over evaluable stations" begin
    fixture = dropout_fixture()
    # Knock the *evaluable* station of fold 1 out at hour 6, leaving the non-evaluable one.
    # Hour 6 then has one evaluable station and it is NaN, so the hour is fully dropped even
    # though the other station's cell is untouched.
    fixture.predictions["mgwr"][2, 6] = NaN
    table = BD.dropout_table(
        fixture.y_obs, fixture.y_sat, fixture.predictions, fixture.ids, fixture.fold_map;
        scheme="balanced_spatial", product="GSMaP", methods=["mgwr"],
    )
    fold1 = only(filter(row -> row.fold == 1, table))
    @test fold1.hours_all_dropped == 2
    @test fold1.cells_in_all_dropped_hours == 3
    @test fold1.cells_scattered == 1
end

@testset "rebuild_common_mask mirrors the benchmark's rule" begin
    y_obs = fill(1.0, 3, 4)
    y_obs[1, 1] = NaN
    a = fill(0.0, 3, 4)
    b = fill(0.0, 3, 4)
    b[2, 2] = NaN
    predictions = Dict{String,Matrix{Float64}}("a" => a, "b" => b, "c" => fill(NaN, 3, 4))
    mask = BD.rebuild_common_mask(y_obs, predictions; methods=["a", "b"])
    @test count(mask) == 10
    @test !mask[1, 1] && !mask[2, 2]
    # A method outside the defining set cannot shrink the mask, however broken it is.
    @test BD.rebuild_common_mask(y_obs, predictions; methods=["a"]) == (.!isnan.(y_obs))
    @test_throws ArgumentError BD.rebuild_common_mask(y_obs, predictions; methods=["missing"])
end

@testset "run_comparison_table pairs the two runs on cells both evaluated" begin
    y_obs = fill(4.0, 3, 4)
    # The "after" run has the wider mask: one extra cell the "before" run had excluded.
    before_mask = trues(3, 4); before_mask[1, 1] = false
    after_mask = trues(3, 4)
    # `steady` predicts identically in both runs, so its paired delta must be exactly zero even
    # though its own-mask RMSE differs — that difference is the mask moving, not the method.
    steady_before = fill(4.0, 3, 4); steady_before[1, 1] = 1.0
    steady_after = copy(steady_before)
    # `improved` halves its error between runs.
    improved_before = fill(6.0, 3, 4)
    improved_after = fill(5.0, 3, 4)
    # `patchy` is NaN at one cell in the after run; the pair must drop it on both sides.
    patchy_before = fill(4.0, 3, 4); patchy_before[3, 3] = 10.0
    patchy_after = fill(4.0, 3, 4); patchy_after[3, 3] = NaN

    before = Dict{String,Matrix{Float64}}(
        "steady" => steady_before, "improved" => improved_before,
        "patchy" => patchy_before, "gone" => fill(4.0, 3, 4),
    )
    after = Dict{String,Matrix{Float64}}(
        "steady" => steady_after, "improved" => improved_after,
        "patchy" => patchy_after, "fresh" => fill(4.0, 3, 4),
    )
    table = BD.run_comparison_table(
        y_obs, before, before_mask, after, after_mask;
        scheme="balanced_spatial", product="GPM",
    )
    @test nrow(table) == 5
    @test all(table.mask_cells_shared .== 11)

    steady = only(filter(:method => ==("steady"), table))
    @test steady.delta_paired == 0.0
    @test steady.n_paired == 11
    # The own-mask numbers do differ, which is exactly the confound the pairing removes.
    @test steady.RMSE_before != steady.RMSE_after

    improved = only(filter(:method => ==("improved"), table))
    @test improved.RMSE_paired_before ≈ 2.0
    @test improved.RMSE_paired_after ≈ 1.0
    @test improved.delta_paired ≈ -1.0
    @test improved.relative_paired ≈ -0.5

    patchy = only(filter(:method => ==("patchy"), table))
    @test patchy.n_paired == 10          # the after-run NaN drops that cell from both sides
    @test patchy.RMSE_paired_before == 0.0

    @test only(filter(:method => ==("gone"), table)).present_in == "before"
    @test only(filter(:method => ==("fresh"), table)).present_in == "after"
    @test isnan(only(filter(:method => ==("gone"), table)).delta_paired)
end

@testset "mask_cost_table reports the cells one method's failures remove" begin
    y_obs = fill(2.0, 3, 4)
    good = fill(2.0, 3, 4)
    broken = fill(2.0, 3, 4)
    broken[1, 1] = NaN
    broken[2, 3] = NaN
    predictions = Dict{String,Matrix{Float64}}("good" => good, "mgwr" => broken)
    table = BD.mask_cost_table(
        y_obs, predictions; scheme="balanced_spatial", product="GPM",
        excluded="mgwr", methods=["good", "mgwr"],
    )
    @test nrow(table) == 2
    @test all(table.mask_cells_full .== 10)
    @test all(table.mask_cells_reduced .== 12)
    @test all(table.cells_recovered .== 2)

    good_row = only(filter(:method => ==("good"), table))
    @test good_row.n_full == 10
    @test good_row.n_reduced == 12          # scored on every cell once mgwr stops defining the mask
    mgwr_row = only(filter(:method => ==("mgwr"), table))
    @test mgwr_row.n_reduced == 10          # its own NaNs are still skipped by the scorer

    @test_throws ArgumentError BD.mask_cost_table(
        y_obs, predictions; scheme="balanced_spatial", product="GPM",
        excluded="good", methods=["mgwr"],
    )
end

"""Write a wide `time × station` table, `time` first, as the benchmark writes its OOF files."""
function write_wide(path, times, ids, matrix)
    df = DataFrame(time=times)
    for (index, id) in enumerate(ids)
        df[!, id] = matrix[index, :]
    end
    CSV.write(path, df)
end

@testset "matched_evaluation_mask: one sample for every product" begin
    ids = ["s1", "s2", "s3"]
    times = [DateTime(2022, 6, 1, 9) + Hour(h) for h in 0:3]
    mktempdir() do scheme_dir
        own = Dict(
            "FY4B" => trues(3, 4), "GPM_uncal" => trues(3, 4), "MERGED_OLS" => trues(3, 4),
        )
        own["FY4B"][1, 2] = false          # a product-specific gap removes the cell for all
        own["GPM_uncal"][3, 4] = false
        blend = ones(3, 4)
        blend[2, 1] = NaN                  # a non-mask method's gap counts too
        for (product, mask) in own
            dir = joinpath(scheme_dir, lowercase(product))
            mkpath(dir)
            write_wide(joinpath(dir, "common_evaluation_mask.csv"), times, ids, mask)
            # `auto` exists for every product; the blend only for GPM_uncal, so an absent file is skipped.
            write_wide(joinpath(dir, "oof_auto.csv"), times, ids, ones(3, 4))
            product == "GPM_uncal" && write_wide(joinpath(dir, "oof_blend_agrenv_mgwr.csv"), times, ids, blend)
            # A mask-defining method's file is not read: its NaN is already in the mask.
            write_wide(joinpath(dir, "oof_mgwr.csv"), times, ids, fill(NaN, 3, 4))
        end
        products = sort(collect(keys(own)))
        sample = BD.matched_evaluation_mask(scheme_dir, products, ids, times)
        expected = trues(3, 4)
        expected[1, 2] = expected[3, 4] = expected[2, 1] = false
        @test sample.matched == expected
        @test sample.masks["FY4B"] == own["FY4B"]
        @test BD.matched_evaluation_mask(scheme_dir, ["MERGED_OLS"], ids, times).matched == trues(3, 4)
        # Files on another hour grid are refused rather than read by position.
        @test_throws ErrorException BD.matched_evaluation_mask(scheme_dir, products, ids, times[2:end])
        @test_throws ArgumentError BD.matched_evaluation_mask(scheme_dir, String[], ids, times)

        table = BD.matched_sample_table(sample.masks, sample.matched, [1, 1, 2]; scheme="balanced_spatial")
        all_rows = filter(:fold => ==("all"), table)
        @test all(all_rows.n_matched .== 9)
        @test only(filter(:product => ==("FY4B"), all_rows)).n_own == 11
        @test only(filter(row -> row.product == "FY4B" && row.fold == "2", table)).n_matched == 3
        @test all(all_rows.hours_matched .== 4)
    end
end

@testset "matched_fold_failures: only failed or skipped folds of the products and methods" begin
    status = DataFrame(
        scheme=["balanced_spatial", "balanced_spatial", "balanced_spatial", "random", "balanced_spatial"],
        product=["FY4B", "GSMaP", "GSMaP", "FY4B", "OTHER"], fold=[1, 2, 3, 1, 1],
        method=["mgwr", "blend_mgwr", "mgwr", "mgwr", "mgwr"],
        status=["partial", "skipped", "failed", "failed", "failed"])
    failures = BD.matched_fold_failures(status; scheme="balanced_spatial", products=["FY4B", "GSMaP"])
    @test failures.fold == [2, 3]
    @test nrow(BD.matched_fold_failures(status; scheme="balanced_spatial", products=["FY4B"])) == 0
end

@testset "--fy4b-file: parsed, and held to the run name's series" begin
    mktempdir() do root
        mkpath(joinpath(root, "output"))
        touch(joinpath(root, "output", "fy4b_w3_h0.csv"))
        path, rest = BD.fy4b_file_option(["run", "--fy4b-file", "output/fy4b_w3_h0.csv", "--x"]; root)
        @test path == joinpath(root, "output", "fy4b_w3_h0.csv") && rest == ["run", "--x"]
        @test BD.fy4b_file_option(["run"]; root) == (nothing, ["run"])
        @test_throws ArgumentError BD.fy4b_file_option(["--fy4b-file"]; root)
        @test_throws ArgumentError BD.fy4b_file_option(["--fy4b-file", "missing.csv"]; root)
    end
    @test BD.check_fy4b_tag("x_fusedanchor_fy4bw3h0_fy4bqm_gpmuncal", "/a/fy4b_w3_h0.csv") === nothing
    @test BD.check_fy4b_tag("x_fusedanchor_nofy4b_gpmuncal", nothing) === nothing
    @test_throws ErrorException BD.check_fy4b_tag("x_fusedanchor_fy4bw3h0_gpmuncal", nothing)
    @test_throws ErrorException BD.check_fy4b_tag("x_fusedanchor_fy4bw3h0_gpmuncal", "/a/fy4b_w1_h0.csv")
    @test_throws ErrorException BD.check_fy4b_tag("x_fusedanchor_gpmuncal", "/a/fy4b_w3_h0.csv")
end
