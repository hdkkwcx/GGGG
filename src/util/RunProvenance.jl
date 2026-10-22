"""
Which code, environment and inputs produced a run, as `(key, value)` rows.

A run directory under `output/` is only citable if it can say what made it. `git_commit` fixes
the code only when `git_dirty` is false, the command line picks what that code does, and `data/`
is outside git so the commit does not fix the inputs. These rows record the git state at the start
(and again at the end, `git_dirty_at_end`), the command, the Julia version, the SHA-256 of the
manifest (copied beside the outputs), and the SHA-256 of every input file.

Plumbing with no opinion about precipitation, shared by the benchmark and the scripts behind the
paper, like `TableIO`. Rows are NamedTuples, so callers write them with `CSV.write` or append
them to a DataFrame without this module needing DataFrames.
"""
module RunProvenance

using SHA: sha256

export start_provenance, finish_provenance!, git_provenance, process_command

const Row = NamedTuple{(:key, :value),Tuple{String,String}}

row(key::AbstractString, value::AbstractString) = Row((String(key), String(value)))

"""
`git_commit`, `git_branch` and `git_dirty` for the checkout at `root`.

Every call is wrapped. A missing git, a checkout that is not a repository, or an ownership check
that refuses records `"unavailable"` rather than aborting a run that takes hours.
"""
function git_provenance(root::AbstractString)
    run_git(args::Vector{String}) = try
        String(strip(read(pipeline(`git -C $root $args`; stderr=devnull), String)))
    catch
        nothing
    end
    porcelain = run_git(["status", "--porcelain"])
    return [
        row("git_commit", something(run_git(["rev-parse", "HEAD"]), "unavailable")),
        row("git_branch", something(run_git(["rev-parse", "--abbrev-ref", "HEAD"]), "unavailable")),
        row("git_dirty", porcelain === nothing ? "unavailable" : string(!isempty(porcelain))),
    ]
end

file_sha256(path::AbstractString) = bytes2hex(open(sha256, path))

"""
The command that started this process: the script's path, relative to `root` when it lies under
it, then `ARGS`, joined with spaces - what follows `julia --project=.` on the command line.
`"unavailable"` when no script is running (the REPL, `julia -e`).
"""
function process_command(root::AbstractString)
    isempty(PROGRAM_FILE) && return "unavailable"
    script = abspath(PROGRAM_FILE)
    root_prefix = rstrip(normpath(root), '/') * "/"
    startswith(script, root_prefix) && (script = script[length(root_prefix)+1:end])
    return join([script; ARGS], ' ')
end

"""
The rows a run records before it does any work: `git_provenance(root)`, `julia_version`,
`manifest_sha256`, `command` (by default `process_command(root)`), and one `input_sha256:<path>`
row per file in `inputs`, the path relative to `root` when it lies under it. The active
environment's `Manifest.toml` is copied into `outdir`.

Read at the start because the code that runs is the code checked out now: read at the end, a run
launched on uncommitted edits that were committed while it ran recorded itself as clean. A dirty
tree is recorded and warned about, never refused. A missing manifest records `"unavailable"`, a
missing input `"missing"`.
"""
function start_provenance(root::AbstractString, outdir::AbstractString; inputs=String[],
                          command::AbstractString=process_command(root))
    rows = git_provenance(root)
    rows[3].value == "true" && @warn "The working tree has uncommitted changes: this run will " *
        "record git_dirty=true, and no commit reproduces it. Commit first for a citable run." root
    push!(rows, row("julia_version", string(VERSION)))
    project = Base.active_project()
    manifest = project === nothing ? nothing : joinpath(dirname(project), "Manifest.toml")
    manifest_sha = if manifest !== nothing && isfile(manifest)
        mkpath(outdir)
        cp(manifest, joinpath(outdir, "Manifest.toml"); force=true)
        file_sha256(manifest)
    else
        "unavailable"
    end
    push!(rows, row("manifest_sha256", manifest_sha))
    push!(rows, row("command", command))
    root_prefix = rstrip(normpath(root), '/') * "/"
    for path in unique(normpath.(inputs))
        name = startswith(path, root_prefix) ? path[length(root_prefix)+1:end] : path
        push!(rows, row("input_sha256:$name", isfile(path) ? file_sha256(path) : "missing"))
    end
    return rows
end

"""Append `git_dirty_at_end`: the second reading, which catches a tree edited during the run."""
function finish_provenance!(rows::AbstractVector, root::AbstractString)
    push!(rows, row("git_dirty_at_end", git_provenance(root)[3].value))
    return rows
end

end # module
