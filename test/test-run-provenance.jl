using Test

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("RunProvenance")
const RP = Main.RunProvenance

@testset "start_provenance records git, Julia, the manifest and every input" begin
    root = mktempdir()
    input = joinpath(root, "data", "gauges.csv")
    mkpath(dirname(input))
    write(input, "time,s1\n2022-01-01T09:00:00,0.5\n")
    outdir = joinpath(root, "output", "run")
    # `root` is not a repository, so the git rows fall back rather than abort.
    rows = RP.start_provenance(root, outdir; inputs=[input, input, joinpath(root, "absent.csv")])
    value(key) = only(r.value for r in rows if r.key == key)
    @test [r.key for r in rows[1:5]] == ["git_commit", "git_branch", "git_dirty", "julia_version", "manifest_sha256"]
    @test value("git_commit") == "unavailable"
    @test value("git_dirty") == "unavailable"
    @test value("julia_version") == string(VERSION)
    @test value("command") == RP.process_command(root)
    # A repeated input is hashed once, under its path relative to `root`.
    @test value("input_sha256:data/gauges.csv") == RP.file_sha256(input)
    @test count(r -> startswith(r.key, "input_sha256:"), rows) == 2
    @test value("input_sha256:absent.csv") == "missing"
    if value("manifest_sha256") != "unavailable"
        @test RP.file_sha256(joinpath(outdir, "Manifest.toml")) == value("manifest_sha256")
    end

    RP.finish_provenance!(rows, root)
    @test rows[end].key == "git_dirty_at_end"

    # An explicit command is recorded as given.
    given = RP.start_provenance(root, outdir; command="scripts/run_x.jl full --flag")
    @test only(r.value for r in given if r.key == "command") == "scripts/run_x.jl full --flag"
end

@testset "process_command is the script, relative to the root, and its arguments" begin
    # Run as a script (this file, or runtests.jl under Pkg.test) the path is under test/; under
    # `julia -e` there is no script.
    command = RP.process_command(joinpath(@__DIR__, ".."))
    @test command == "unavailable" || startswith(command, "test/")
end

@testset "git_provenance reads this checkout" begin
    rows = RP.git_provenance(joinpath(@__DIR__, ".."))
    @test [r.key for r in rows] == ["git_commit", "git_branch", "git_dirty"]
    # In a checkout without git available every row is "unavailable"; otherwise a 40-hex commit.
    commit = rows[1].value
    @test commit == "unavailable" || occursin(r"^[0-9a-f]{40}$", commit)
    @test rows[3].value in ("true", "false", "unavailable")
end
