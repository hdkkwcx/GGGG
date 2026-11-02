using Test, DataFrames, Dates, NCDatasets

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("FY4BCloudType")
const CLT = Main.FY4BCloudType
const FY4B_PRE = Main.FY4BPreprocessing

const CLT_STATION = (lat=32.3, lon=110.5)

clt_name(start; product="CLT") = "FY4B-_AGRI--_N_DISK_1330E_L2-_$(product)-_MULT_NOM_" *
    Dates.format(start, "yyyymmddHHMMSS") * "_" *
    Dates.format(start + Minute(15) - Second(1), "yyyymmddHHMMSS") * "_4000M_V0001.NC"

"""A synthetic full-disk CLT file: all water (2), with `paint` applied as (x, y) => code."""
function write_clt(path; paint=Dict{Tuple{Int,Int},Int}())
    codes = fill(Int8(2), FY4B_PRE.GRID_SIZE, FY4B_PRE.GRID_SIZE)
    for ((x, y), code) in paint
        codes[x, y] = Int8(code)
    end
    NCDataset(path, "c") do ds
        defDim(ds, "y", FY4B_PRE.GRID_SIZE); defDim(ds, "x", FY4B_PRE.GRID_SIZE)
        defVar(ds, "CLT", codes, ("y", "x"); deflatelevel=1)
    end
end

@testset "CLT files: strict hours, Beijing labels, QPE names ignored" begin
    mktempdir() do dir
        hour1 = DateTime(2023, 7, 1, 0)           # complete
        hour2 = DateTime(2023, 7, 1, 1)           # three quarters only
        for minute_ in (0, 15, 30, 45)
            write_clt(joinpath(dir, clt_name(hour1 + Minute(minute_))))
        end
        for minute_ in (0, 15, 45)
            write_clt(joinpath(dir, clt_name(hour2 + Minute(minute_))))
        end
        touch(joinpath(dir, clt_name(hour2 + Minute(30); product="QPE")))  # not a CLT file
        touch(joinpath(dir, "._" * clt_name(hour2 + Minute(30))))           # AppleDouble companion
        groups, qc = CLT.find_clt_files(dir; years=2023:2023, months=7:7)
        @test length(groups[hour1]) == 4 && length(groups[hour2]) == 3
        @test nrow(qc) == 31 * 24
        row1 = only(eachrow(qc[qc.hour_start .== hour1, :]))
        row2 = only(eachrow(qc[qc.hour_start .== hour2, :]))
        @test row1.complete && row1.time == DateTime(2023, 7, 1, 9)        # UTC 00-01 -> 09:00 BJT
        @test !row2.complete && row2.missing_minutes == "30"
        stations = DataFrame(lat=[CLT_STATION.lat], lon=[CLT_STATION.lon])
        result = CLT.clt_class_counts(groups, stations, [(0, 0)])
        @test result.hour_start == [hour1] && result.time == [DateTime(2023, 7, 1, 9)]
        water = findfirst(==(2), CLT.CLT_CLASSES)
        @test result.counts[1, 1, 1, water] == 36 && sum(result.counts) == 36
    end
end

@testset "CLT box counts follow the shift; space and fill are not classes" begin
    mktempdir() do dir
        x, y = FY4B_PRE._pixel(FY4B_PRE.latlon_to_xy(CLT_STATION.lat, CLT_STATION.lon, 133.0)...)
        hour = DateTime(2023, 7, 2, 6)
        # Every quarter: cirrus (6) on the nominal centre; ice (5) one pixel east and one south of
        # it, then space (126) on the pixel east of that.
        paint = Dict((x, y) => 6, (x + 1, y + 1) => 5, (x + 2, y + 1) => 126)
        for minute_ in (0, 15, 30, 45)
            write_clt(joinpath(dir, clt_name(hour + Minute(minute_))); paint)
        end
        # The last quarter also has fill (127) at the centre.
        rm(joinpath(dir, clt_name(hour + Minute(45))))
        write_clt(joinpath(dir, clt_name(hour + Minute(45)));
            paint=merge(paint, Dict((x, y) => 127)))
        groups, _ = CLT.find_clt_files(dir; years=2023:2023, months=7:7)
        stations = DataFrame(lat=[CLT_STATION.lat], lon=[CLT_STATION.lon])
        shifts = [(0, 0), (1, 1)]
        counts = CLT.clt_class_counts(groups, stations, shifts).counts
        index(code) = findfirst(==(code), CLT.CLT_CLASSES)
        # Unshifted box: the centre is cirrus in 3 quarters (fill in the 4th); ice at (+1, +1).
        @test counts[1, 1, 1, index(6)] == 3 && counts[1, 1, 1, index(5)] == 4
        @test counts[1, 1, 1, index(2)] == 4 * 7 && sum(counts[1, 1, :, :]) == 35
        # Shifted (+1, +1): the box is centred on the ice pixel; (x+2, y+1) is space.
        @test counts[1, 2, 1, index(5)] == 4 && counts[1, 2, 1, index(6)] == 3
        @test sum(counts[1, 2, :, :]) == 36 - 4 - 1
        # Dominant: water wins both boxes; with too few valid pixel-quarters, -1.
        dominant = CLT.dominant_class(counts; min_valid=18)
        @test dominant[1, 1, 1] == 2 && dominant[1, 2, 1] == 2
        @test all(==(-1), CLT.dominant_class(counts; min_valid=36))
    end
end

@testset "dominant class: a tie goes to the smaller code" begin
    counts = zeros(UInt8, 1, 1, 2, length(CLT.CLT_CLASSES))
    index(code) = findfirst(==(code), CLT.CLT_CLASSES)
    counts[1, 1, 1, index(6)] = 18; counts[1, 1, 1, index(5)] = 18     # cirrus ties ice
    counts[1, 1, 2, index(7)] = 10; counts[1, 1, 2, index(0)] = 3      # overlap, too few valid
    dominant = CLT.dominant_class(counts; min_valid=18)
    @test dominant[1, 1, 1] == 5 && dominant[1, 1, 2] == -1
end
