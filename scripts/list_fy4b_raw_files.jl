#!/usr/bin/env julia

# List the raw FY4B 15-minute QPE files on the external drive: where each copy lies, how large it
# is, and which copy the rebuild's rule keeps. The raw tree is the one input of the paper's run with
# no SHA-256, and the rebuild's fy4b_hourly_qc.csv keeps file names only, so this is the record of
# what the names stood for. It checks names and sizes, not every byte: nothing is hashed.
#
#   julia --project=. scripts/list_fy4b_raw_files.jl     # -> output/fy4b_raw_manifest/
#
# The filters are the rebuild's own (`FY4BPreprocessing`): AppleDouble companions and files that
# are not 15-minute scans are skipped, so the gauge-corrected 1-, 3- and 24-hour files never appear.
# Of several copies of one scan, `find_nc_files` keeps the first that opens, in sorted path order.
# Opening 29,007 repeated copies again would take the rebuild's time, so here only the scans whose
# copies differ in size are opened (a damaged download is the usual reason): `kept` is the first
# copy that opens there, and the first in sorted path order everywhere else. The kept names are
# checked against the QC table.
# The drive is only read. scripts/trace_fy4b_process.jl carries summary.csv and checks.csv into the
# FY4B process page when they exist.

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("FY4BPreprocessing", "RunProvenance", "TableIO")

using CSV, DataFrames, Dates
using Main.FY4BPreprocessing: is_appledouble, is_15min_file, parse_nc_filename, segment_problem
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic

const RAW_DIR = joinpath(ROOT, "data", "FY4B")
const QC_FILE = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024", "fy4b_hourly_qc.csv")
const REBUILD_LOG = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024.log")
const OUTDIR = joinpath(ROOT, "output", "fy4b_raw_manifest")
const YEARS = 2022:2024

function list_files()
    rows = NamedTuple[]
    for (dir, _, files) in walkdir(RAW_DIR; follow_symlinks=true), file in files
        is_appledouble(file) && continue
        endswith(uppercase(file), ".NC") || continue
        is_15min_file(file) || continue
        start, _, sat_lon = parse_nc_filename(file)
        year(start) in YEARS || continue
        path = joinpath(dir, file)
        relative = relpath(path, RAW_DIR)
        info = stat(path)
        push!(rows, (; folder=first(splitpath(relative)), path=relative, name=file, scan_start=start, sat_lon,
            bytes=Int(info.size), modified=Dates.format(unix2datetime(info.mtime), "yyyy-mm-ddTHH:MM:SS")))
    end
    table = DataFrame(rows)
    sort!(table, [:scan_start, :sat_lon, :path])
    table.copies = zeros(Int, nrow(table))
    table.kept = falses(nrow(table))
    table.same_size_as_kept = trues(nrow(table))
    # "" where the copy was not opened: every copy of the scan has one size.
    table.opens = fill("", nrow(table))
    for group in groupby(table, [:scan_start, :sat_lon])
        group.copies .= nrow(group)
        keep = 1
        if !allequal(group.bytes)
            opens = [segment_problem(joinpath(RAW_DIR, path)) === nothing for path in group.path]
            group.opens .= string.(opens)
            keep = something(findfirst(opens), 1)
        end
        group.kept[keep] = true
        group.same_size_as_kept .= group.bytes .== group.bytes[keep]
    end
    return table
end

function main()
    isdir(RAW_DIR) || error("Raw FY4B tree not found at $RAW_DIR - is the external drive mounted?")
    provenance = start_provenance(ROOT, OUTDIR; inputs=[QC_FILE])
    table = list_files()
    kept = table[table.kept, :]
    repeated = table[table.copies .> 1, :]

    qc = CSV.read(QC_FILE, DataFrame; types=Dict(:source_files => String, :unreadable_files => String))
    used = String[]
    for names in skipmissing(qc.source_files), name in split(names, ';'; keepempty=false)
        push!(used, name)
    end
    unreadable = Set(basename(name) for names in skipmissing(qc.unreadable_files)
                     for name in split(names, ';'; keepempty=false))
    table.unreadable_in_rebuild = in(unreadable).(table.name)

    folders = combine(groupby(table, :folder), nrow => :files, :kept => sum => :kept_files,
        :bytes => sum => :bytes, :scan_start => minimum => :first_scan, :scan_start => maximum => :last_scan)
    sort!(folders, :folder)
    summary = DataFrame(item=String[], value=String[])
    push!(summary, ("15-minute files on the drive, every copy", string(nrow(table))))
    push!(summary, ("distinct scans (start time and satellite position)", string(nrow(kept))))
    push!(summary, ("repeated copies", string(nrow(table) - nrow(kept))))
    push!(summary, ("scans with more than one copy", string(count(>(1), kept.copies))))
    push!(summary, ("scans whose copies differ in size", string(count(
        group -> !allequal(group.bytes), groupby(table, [:scan_start, :sat_lon])))))
    push!(summary, ("copies that differ in size from the kept copy", string(count(.!table.same_size_as_kept))))
    push!(summary, ("of them, copies that do not open (damaged downloads)",
        string(count(.!table.same_size_as_kept .& (table.opens .== "false")))))
    push!(summary, ("size of the kept files, GB", string(round(sum(kept.bytes) / 1e9; digits=1))))
    push!(summary, ("size of every copy, GB", string(round(sum(table.bytes) / 1e9; digits=1))))
    push!(summary, ("smallest and largest kept file, MB",
        "$(round(minimum(kept.bytes) / 1e6; digits=2)) - $(round(maximum(kept.bytes) / 1e6; digits=2))"))
    for row in eachrow(folders)
        push!(summary, ("folder $(row.folder): files, of which kept", "$(row.files), $(row.kept_files)"))
    end

    checks = DataFrame(check=String[], expected=String[], found=String[], pass=Bool[])
    check!(name, expected, found) = push!(checks, (name, string(expected), string(found), expected == found))
    check!("raw drive: distinct scans = files the rebuild used (QC table)", sum(qc.file_count), nrow(kept))
    check!("raw drive: the kept copies' names = the QC table's source files, name for name",
        sort(used), sort(kept.name))
    if isfile(REBUILD_LOG)
        found = match(r"duplicate_file_count = (\d+)", read(REBUILD_LOG, String))
        found === nothing || check!("raw drive: repeated copies = the rebuild log's count",
            parse(Int, found.captures[1]), nrow(table) - nrow(kept))
    end
    check!("raw drive: every copy that differs in size from the kept copy is a file that does not open", 0,
        count(.!table.same_size_as_kept .& (table.opens .!= "false")))
    check!("raw drive: where copies differ in size, the kept copy opens", 0,
        count(table.kept .& (table.opens .== "false")))
    check!("raw drive: no kept file is empty", 0, count(==(0), kept.bytes))
    # The check on names compares long vectors; keep the table readable.
    checks.expected = [length(text) > 200 ? "$(length(used)) names" : text for text in checks.expected]
    checks.found = [length(text) > 200 ? "$(nrow(kept)) names" : text for text in checks.found]

    write_csv_atomic(joinpath(OUTDIR, "raw_files.csv"), table)
    write_csv_atomic(joinpath(OUTDIR, "copies.csv"), repeated)
    write_csv_atomic(joinpath(OUTDIR, "folders.csv"), folders)
    write_csv_atomic(joinpath(OUTDIR, "summary.csv"), summary)
    write_csv_atomic(joinpath(OUTDIR, "checks.csv"), checks)
    CSV.write(joinpath(OUTDIR, "provenance.csv"), finish_provenance!(provenance, ROOT))
    show(summary; allrows=true, truncate=0)
    println()
    show(checks; allrows=true, truncate=0)
    println()
    all(checks.pass) || error("$(count(.!checks.pass)) raw-file check(s) failed; see checks.csv")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
