#!/usr/bin/env julia

# Rebuild the paper's figures and cited numbers in dependency order.
#
#   julia --project=. scripts/make_paper.jl                 # dry run: what is reused, what would rerun
#   julia --project=. scripts/make_paper.jl --only fig06b   # one figure, table or step, and what it needs
#   julia --project=. scripts/make_paper.jl --run           # rerun what is missing or stale, refresh paper/
#
# Walks the `[[step]]` graph in paper/sources.toml. An analysis or figure step is REUSE while its
# marker exists and, when the marker is a provenance file, records a clean tree (`git_dirty*` =
# false, as scripts/check_paper_provenance.py requires) and inputs whose SHA-256 still match; STALE
# when it does not, or when a step it needs will rerun; MISSING without a marker. Data steps write
# data/, so they are only checked and their command printed, never run: OK, MISSING, or ABSENT for
# an intermediate that was cleared after the data steps needing it were made.
#
# A REUSE step is CHANGED when the code it ran has changed since the commit its provenance records:
# for a Julia step, anything under src/, Project.toml, its script or the Manifest; for a plot step,
# its script, scripts/plot_provenance.py or scripts/requirements-plot.txt. CHANGED is reported, never
# rerun on its own, since most edits do not move a result and the benchmark takes hours: review the
# diff, then rerun with `--force`, or record the review as `reviewed = "<commit>"` on the step, after
# which only changes past that commit count. A step is also CHANGED when its provenance records a
# `command` (runs since 2026-10-21 do) other than the step's own, from its `scripts/` path on.
#
# `--run` refuses on a dirty tree, as every paper run starts clean. It runs the STALE and MISSING
# steps with their own commands (logs in output/logs/<step>.log), copies each figure into
# paper/figures/ when it differs from its source, re-exports paper/numbers/, runs the provenance
# check and compiles paper/paper_plan.typ when typst is installed. A dry run only runs the check.
#
# `--force <step>` treats one step as STALE; `--sources <path>` reads another sources.toml (for
# testing). Python is $PYTHON, else python3.13, else python3: the paper scripts need Python 3.11+
# (tomllib), and the plot scripts also need scripts/requirements-plot.txt installed.

const ROOT = normpath(joinpath(@__DIR__, ".."))

using SHA, TOML

const CLEAN_KEYS = ("git_dirty", "git_dirty_at_end", "report_git_dirty")

file_sha256(path) = open(io -> bytes2hex(sha256(io)), path)

"""The `key => value` rows of a provenance CSV (header `key,value`)."""
read_provenance(path) = Dict(String(first(row)) => String(last(row))
    for row in (split(line, ','; limit=2) for line in Iterators.drop(eachline(path), 1))
    if length(row) == 2)

function python()
    haskey(ENV, "PYTHON") && return split(ENV["PYTHON"])
    for name in ("python3.13", "python3")
        path = Sys.which(name)
        path === nothing || return [path]
    end
    error("no Python found; set PYTHON")
end

"""A step's command, with "julia" as the running Julia and "python" as `python()`."""
function step_cmd(step)
    argv = step["command"]
    head = argv[1] == "julia" ? Base.julia_cmd().exec : argv[1] == "python" ? python() : [argv[1]]
    return Cmd(Cmd(String[head; argv[2:end]]); dir=ROOT)
end

"""`(status, reason)` of one step from its own marker or outputs, before looking at its needs."""
function marker_status(step)
    if step["kind"] == "data"
        absent = filter(path -> !ispath(joinpath(ROOT, path)), step["outputs"])
        return isempty(absent) ? ("OK", "") : ("MISSING", "absent: " * join(absent, ", "))
    end
    marker = joinpath(ROOT, step["marker"])
    isfile(marker) || return ("MISSING", "no $(step["marker"])")
    readline(marker) == "key,value" || return ("REUSE", "code not checked: no commit recorded")
    record = read_provenance(marker)
    dirty = [key for key in CLEAN_KEYS if get(record, key, "false") != "false"]
    isempty(dirty) || return ("STALE", join(dirty, ", ") * " in its provenance")
    for (key, recorded) in record
        startswith(key, "input_sha256:") || continue
        relative = chopprefix(key, "input_sha256:")
        path = isabspath(relative) ? relative : joinpath(ROOT, relative)
        isfile(path) || return ("STALE", "input gone: $relative")
        file_sha256(path) == recorded || return ("STALE", "input changed: $relative")
    end
    # Runs since 2026-10-21 record their command; the step's must be the one that made the result.
    recorded = get(record, "report_command", get(record, "command", nothing))
    argv = step["command"]
    expected = join(argv[something(findfirst(startswith("scripts/"), argv), 1):end], ' ')
    recorded === nothing || recorded == expected ||
        return ("CHANGED", "made by a different command: $recorded")
    base, changed, note = code_changes(step, record)
    isempty(changed) && return ("REUSE", note)
    return ("CHANGED", strip("code changed since $(base[1:7]): $(join(changed, ", ")). $note"))
end

git_ok(args...) = success(pipeline(Cmd(`git $args`; dir=ROOT); stderr=devnull))

"""
`(base, changed, note)`: the files of the step's code that changed between `base` and the
working tree, plus any environment that no longer matches the recorded one. `base` is the commit
the step ran at, or its `reviewed` commit when that descends from it: code changes up to `reviewed`
were checked and do not move the result. Julia steps count `src/` as a whole (the loaders pull
modules in by name), `Project.toml`, their script and the Manifest; plot steps their script,
`scripts/plot_provenance.py` and `scripts/requirements-plot.txt`.
"""
function code_changes(step, record)
    commit = get(record, "report_git_commit", get(record, "git_commit", ""))
    isempty(commit) && return ("", String[], "code not checked: no commit recorded")
    git_ok("cat-file", "-e", "$commit^{commit}") ||
        return (commit, ["recorded commit $(commit[1:7]) is not in this repository"], "")
    base, note = commit, ""
    if haskey(step, "reviewed")
        reviewed = step["reviewed"]
        if git_ok("merge-base", "--is-ancestor", commit, reviewed)
            base = reviewed
        elseif !git_ok("merge-base", "--is-ancestor", reviewed, commit)   # older: the run supersedes it
            note = "reviewed $(reviewed[1:min(7, end)]) ignored: not a descendant of $(commit[1:7])"
        end
    end
    julia = step["command"][1] == "julia"
    scripts = filter(startswith("scripts/"), step["command"])
    paths = julia ? ["src", "Project.toml", scripts...] : [scripts..., "scripts/plot_provenance.py"]
    changed = split(read(Cmd(`git diff --name-only $base -- $paths`; dir=ROOT), String))
    environment = julia ? ("manifest_sha256", "Manifest.toml") :
        ("requirements_sha256", "scripts/requirements-plot.txt")
    recorded = get(record, environment[1], nothing)
    path = joinpath(ROOT, environment[2])
    if recorded !== nothing && (!isfile(path) || file_sha256(path) != recorded)
        push!(changed, "$(environment[2]) (differs from the recorded one)")
    end
    return (base, String.(changed), note)
end

"""Step ids in dependency order, data then analysis then figures, failing on an unknown id or a cycle."""
function dependency_order(byid)
    order, state = String[], Dict{String,Symbol}()
    function visit(id, from)
        haskey(byid, id) || error("step $from needs unknown step $id")
        get(state, id, :new) == :done && return
        get(state, id, :new) == :active && error("dependency cycle through step $id")
        state[id] = :active
        foreach(need -> visit(need, id), get(byid[id], "needs", String[]))
        state[id] = :done
        push!(order, id)
    end
    foreach(id -> visit(id, id), sort!(collect(keys(byid))))
    # Stable, and still a dependency order: data needs only data, analysis no figure.
    rank = Dict("data" => 1, "analysis" => 2, "figure" => 3)
    return sort!(order; by=id -> rank[byid[id]["kind"]], alg=MergeSort)
end

"""`id` and every step it needs, directly or not."""
function ancestors(byid, id, found=Set{String}())
    id in found && return found
    push!(found, id)
    foreach(need -> ancestors(byid, need, found), get(byid[id], "needs", String[]))
    return found
end

function parse_args(args)
    options = Dict{String,Any}("run" => false, "only" => nothing, "force" => String[],
        "sources" => joinpath(ROOT, "paper", "sources.toml"))
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--run"
            options["run"] = true
        elseif arg in ("--only", "--force", "--sources") && i < length(args)
            i += 1
            arg == "--force" ? push!(options["force"], args[i]) : (options[arg[3:end]] = args[i])
        else
            error("unknown argument $arg; see the header of scripts/make_paper.jl")
        end
        i += 1
    end
    return options
end

function main(args=ARGS)
    options = parse_args(args)
    sources = TOML.parsefile(options["sources"])
    items, tables = sources["item"], get(sources, "table", Any[])
    byid = Dict(step["id"] => step for step in sources["step"])
    order = dependency_order(byid)
    for id in options["force"]
        haskey(byid, id) || error("--force: no step $id")
    end

    selected = Set(order)
    if options["only"] !== nothing
        target = options["only"]
        entry = findfirst(e -> get(e, "id", get(e, "name", nothing)) == target, [items; tables])
        step_id = entry !== nothing ? [items; tables][entry]["step"] :
            haskey(byid, target) ? target : error("--only: no figure, table or step $target")
        selected = ancestors(byid, step_id)
    end

    status = Dict{String,Tuple{String,String}}()
    for id in filter(in(selected), order)
        step = byid[id]
        status[id] = id in options["force"] ? ("STALE", "forced") : marker_status(step)
        if step["kind"] != "data" && status[id][1] in ("REUSE", "CHANGED")
            rerun = findfirst(need -> byid[need]["kind"] != "data" &&
                status[need][1] in ("STALE", "MISSING"), get(step, "needs", String[]))
            rerun === nothing || (status[id] = ("STALE", "needs $(step["needs"][rerun]), which reruns"))
        end
    end
    # An intermediate a data step cleared after use is no gap while every data step that needs it
    # is present: it would only matter for remaking them.
    for id in reverse(filter(in(selected), order))
        status[id][1] == "MISSING" && byid[id]["kind"] == "data" || continue
        users = filter(other -> id in get(byid[other], "needs", String[]), collect(selected))
        if !isempty(users) && all(user -> byid[user]["kind"] == "data" && status[user][1] == "OK", users)
            status[id] = ("ABSENT", "only needed to remake $(join(sort(users), ", ")), which is present")
        end
    end

    shown = startswith(abspath(options["sources"]), ROOT) ? relpath(options["sources"], ROOT) : options["sources"]
    println("Steps in dependency order ($shown):")
    to_run = String[]
    for id in filter(in(selected), order)
        step, (state, reason) = byid[id], status[id]
        println("  ", rpad(state, 8), rpad(id, 20), rpad(step["kind"], 9), reason)
        if state in ("STALE", "MISSING")
            command = isempty(step["command"]) ? "(no tracked command)" : join(step["command"], ' ')
            println(" "^10, "command: ", command)
            haskey(step, "note") && println(" "^10, "note: ", step["note"])
            step["kind"] == "data" || push!(to_run, id)
        end
    end
    changed = filter(id -> status[id][1] == "CHANGED", filter(in(selected), order))
    isempty(changed) || println("""
        $(length(changed)) step(s) reused although code they ran has changed since: $(join(changed, ", ")).
        Not rerun automatically. Review each diff (git diff <commit> -- <files>); if a change can move
        the result, rerun with --run --force <step>, otherwise record the review as `reviewed = "<commit>"`
        on the step in paper/sources.toml.""")

    if options["run"]
        dirty = read(Cmd(`git status --porcelain`; dir=ROOT), String)
        isempty(dirty) || error("--run needs a clean tree, as every paper run does; commit or stash:\n$dirty")
        for id in to_run
            blocked = filter(need -> status[need][1] == "MISSING" && byid[need]["kind"] == "data",
                get(byid[id], "needs", String[]))
            isempty(blocked) || error("step $id needs data steps that are missing: $(join(blocked, ", "))")
            log = joinpath(ROOT, "output", "logs", id * ".log")
            mkpath(dirname(log))
            println("Running $id (log: $(relpath(log, ROOT)))")
            flush(stdout)
            open(io -> run(pipeline(step_cmd(byid[id]); stdout=io, stderr=io)), log, "w")
            state, reason = marker_status(byid[id])
            state == "REUSE" || error("step $id ran, but its marker is $state: $reason")
            status[id] = (state, reason)
        end
    else
        println(isempty(to_run) ? "Nothing to rerun." :
            "Would rerun: $(join(to_run, ", ")) (pass --run).")
    end

    copies = 0
    for item in filter(item -> item["step"] in selected, items)
        source, paper = joinpath(ROOT, item["source"]), joinpath(ROOT, item["paper"])
        isfile(source) || (println("  no source for $(item["id"]): $(item["source"])"); continue)
        isfile(paper) && file_sha256(paper) == file_sha256(source) && continue
        copies += 1
        options["run"] ? cp(source, paper; force=true) : println("  would copy $(item["source"]) -> $(item["paper"])")
    end
    println(options["run"] ? "Copied $copies figure(s) into paper/figures/." :
        "$copies figure(s) differ from their source.")

    py = python()
    flush(stdout)
    options["run"] && run(Cmd(Cmd([py; "scripts/export_paper_numbers.py"]); dir=ROOT))
    ok = success(run(ignorestatus(Cmd(Cmd([py; "scripts/check_paper_provenance.py"]); dir=ROOT))))
    if options["run"] && ok
        typst = Sys.which("typst")
        typst === nothing ? println("typst not found; paper/paper_plan.typ not compiled") :
            run(Cmd(`$typst compile paper/paper_plan.typ`; dir=ROOT))
    end
    return ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main() ? 0 : 1)
end
