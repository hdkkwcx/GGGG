# ---------------------------------------------------------------------------- reading

"""
Where a diagnosis of `run_dir` is written: `output/benchmark_diagnostics/<run name>`.

Each diagnostic run gets its own subdirectory named for the benchmark run it read, so diagnosing a
second run never overwrites the first one's numbers. Repeated verbatim in
`run_benchmark_diagnostics.jl`, `run_mgwr_diagnostics.jl` and `run_claim_reassessment.jl`, which
also each checked `run_dir` exists first - done here so the error is worded once.
"""
function benchmark_diagnostics_outdir(root::AbstractString, run_dir::AbstractString)
    isdir(run_dir) || error("benchmark output directory not found: $run_dir")
    outdir = joinpath(root, "output", "benchmark_diagnostics", basename(run_dir))
    mkpath(outdir)
    return outdir
end

"""
`(fy4b_file, rest)`: a script's `--fy4b-file <path>` option, the FY4B series a benchmark run was
given with the same flag (a relative path is read from `root`), and the arguments without it.
"""
function fy4b_file_option(args; root::AbstractString)
    rest = collect(String, args)
    index = findfirst(==("--fy4b-file"), rest)
    index === nothing && return nothing, rest
    index < length(rest) || throw(ArgumentError("--fy4b-file needs a path"))
    path = isabspath(rest[index + 1]) ? rest[index + 1] : joinpath(root, rest[index + 1])
    isfile(path) || throw(ArgumentError("--fy4b-file: no such file $path"))
    return path, deleteat!(rest, index:(index + 1))
end

"""
Refuse to read a run on another FY4B series' grid: the run name carries a rebuilt series' tag
(`fy4b_w3_h0.csv` -> `_fy4bw3h0`, as `run_interpolation_benchmark.jl` names it) exactly when
`fy4b_file` is that series.
"""
function check_fy4b_tag(run_name::AbstractString, fy4b_file::Union{Nothing,AbstractString})
    tagged = match(r"_fy4bw\d+h\d+(?=_|$)", run_name)
    if fy4b_file === nothing
        tagged === nothing ||
            error("$run_name was run on a rebuilt FY4B series; name it with --fy4b-file")
        return nothing
    end
    tag = "_" * replace(splitext(basename(fy4b_file))[1], "_" => "")
    tagged !== nothing && tagged.match == tag ||
        error("--fy4b-file $(basename(fy4b_file)) does not match the run $run_name (expected $tag)")
    return nothing
end

"""Read a wide `time × station` benchmark CSV into a `station × time` matrix ordered by `ids`."""
function read_wide_matrix(path::AbstractString, ids::Vector{String})
    df = CSV.read(path, DataFrame)
    out = Matrix{Float64}(undef, length(ids), nrow(df))
    for (index, id) in enumerate(ids)
        column = df[!, Symbol(id)]
        out[index, :] = [value === missing ? NaN : Float64(value) for value in column]
    end
    return out
end

"""Read `common_evaluation_mask.csv` into a `station × time` `BitMatrix` ordered by `ids`."""
function read_mask_matrix(path::AbstractString, ids::Vector{String})
    df = CSV.read(path, DataFrame)
    out = falses(length(ids), nrow(df))
    for (index, id) in enumerate(ids)
        out[index, :] = Bool.(df[!, Symbol(id)])
    end
    return out
end

"""Read every `oof_<method>.csv` in one `(scheme, product)` directory into `method => matrix`."""
function load_prediction_matrices(product_dir::AbstractString, ids::Vector{String})
    predictions = Dict{String,Matrix{Float64}}()
    for path in readdir(product_dir; join=true)
        name = basename(path)
        startswith(name, "oof_") && endswith(name, ".csv") || continue
        predictions[name[5:end-4]] = read_wide_matrix(path, ids)
    end
    return predictions
end

"""Read `split_common.csv` into `station_id => fold`."""
function read_fold_map(path::AbstractString)
    df = CSV.read(path, DataFrame; types=Dict(:station_id => String))
    return Dict(String(row.station_id) => Int(row.fold) for row in eachrow(df))
end

