# The formulas of the FY4B methods page (scripts/plot_fy4b_methods.py), each written out plainly
# here and compared with the function it describes, on small synthetic inputs. The page states a
# formula as "checked by number" only if its testset below exists and passes; when the code
# changes, the failing testset says which formula on the page is out of date.
#
# The transcriptions share no code with src/: loops and textbook linear algebra only. Names carry
# an `fm_` prefix because runtests.jl includes every test file into one process.

using Test, DataFrames, Dates, LinearAlgebra, Random, Statistics
using MixedGWR

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_pipeline("InterpolationBenchmark")
load_standalone_modules("FY4BPreprocessing", "SatelliteCalibration", "TraditionalInterpolation",
    "JointCovariateModels", "JointVariableSelection", "SatelliteFusion", "DEMTerrainExperiment")

const FM_FY4B = Main.FY4BPreprocessing
const FM_SC = Main.SatelliteCalibration
const FM_TI = Main.TraditionalInterpolation
const FM_JCM = Main.JointCovariateModels
const FM_JVS = Main.JointVariableSelection
const FM_SF = Main.SatelliteFusion
const FM_DEM = Main.DEMTerrainExperiment
const FM_R_KM = 6371.0088

fm_close(a, b; tol=1e-10) = size(a) == size(b) &&
    all((isnan(x) && isnan(y)) || isapprox(x, y; atol=tol, rtol=tol) for (x, y) in zip(a, b))

# --- A. Raw file to hourly value ---------------------------------------------------------------

"""Column and row (0-based, fractional) of a point, as the page writes the navigation."""
function fm_navigation(lat, lon, sat_lon; h=0.0)
    r_e, r_p, H, N, res = 6378.137e3, 6356.75231414e3, 42164.0e3, 2748, 4000.0
    phi, lam = lat * pi / 180, (lon - sat_lon) * pi / 180
    f = (r_e - r_p) / r_e
    e2 = f * (2 - f)
    phi_c = atan((1 - e2) * tan(phi))
    r_c = r_p / sqrt(1 - e2 * cos(phi_c)^2)
    x = r_c * cos(phi_c) * cos(lam) + h * cos(phi) * cos(lam)
    y = r_c * cos(phi_c) * sin(lam) + h * cos(phi) * sin(lam)
    z = r_c * sin(phi_c) + h * sin(phi)
    alpha_x = atan(y / (H - x))
    alpha_y = atan(z / sqrt((H - x)^2 + y^2))
    delta = res / (H - r_e)
    return (N - 1) / 2 + alpha_x / delta, (N - 1) / 2 - alpha_y / delta
end

@testset "formula: navigation" begin
    for (lat, lon) in ((32.3, 110.5), (31.25, 110.312), (33.1, 109.4)), sat_lon in (133.0, 105.0)
        @test all(isapprox.(fm_navigation(lat, lon, sat_lon), FM_FY4B.latlon_to_xy(lat, lon, sat_lon); atol=1e-9))
        @test all(isapprox.(fm_navigation(lat, lon, sat_lon; h=10_000.0),
            FM_FY4B.latlon_to_xy(lat, lon, sat_lon; height=10_000.0); atol=1e-9))
    end
    # The sub-satellite point is the grid centre.
    @test all(isapprox.(fm_navigation(0.0, 133.0, 133.0), (1373.5, 1373.5); atol=1e-9))
end

"""Mean of the accepted pixels of the 3x3 box, and the share of them at the 20 mm/h ceiling."""
function fm_window(block, quality, ix, iy)
    accepted = Float64[]
    for jy in (iy - 1):(iy + 1), jx in (ix - 1):(ix + 1)
        v, q = block[jx, jy], quality[jx, jy]
        (v != Float32(65534.0) && 0 <= v <= 30 && q < 2) && push!(accepted, Float64(v))
    end
    length(accepted) >= 5 || return (NaN, NaN)
    return sum(accepted) / length(accepted), count(>=(20.0), accepted) / length(accepted)
end

@testset "formula: pixel rule and 3x3 mean" begin
    rng = MersenneTwister(1)
    for trial in 1:200
        block = Float32.(rand(rng, 5, 5) .* 22)
        block[rand(rng, 1:25, 3)] .= rand(rng, Float32[65533, 65534, 65535, 20, 20, 31])
        quality = rand(rng, 0:3, 5, 5) .* (rand(rng, 5, 5) .< 0.3)
        ours, theirs = fm_window(block, quality, 3, 3), FM_FY4B._window_value(block, quality, 3, 3, 3)
        @test fm_close(collect(ours), collect(theirs))
    end
end

# --- D. Corrected FY4B --------------------------------------------------------------------------

"""Pearson r of x and y from n, Σx, Σy, Σx², Σy², Σxy."""
function fm_pooled_r(s)
    n = s[1]
    mx, my = s[2] / n, s[3] / n
    return (s[6] / n - mx * my) / sqrt((s[4] / n - mx^2) * (s[5] / n - my^2))
end

@testset "formula: pooled correlation and the shift" begin
    rng = MersenneTwister(2)
    x, y = rand(rng, 400), rand(rng, 400)
    sums = [400.0, sum(x), sum(y), sum(abs2, x), sum(abs2, y), sum(x .* y)]
    @test fm_pooled_r(sums) ≈ cor(x, y)
    @test FM_SC._pooled_r(sums) ≈ cor(x, y)

    shifts = [(dx, dy) for dy in -1:1 for dx in -1:1]
    grid = FM_SC.OffsetGrid("", shifts, 3, [1, 2, 3], [1])
    table = zeros(length(shifts), 3, 1, 6)
    for k in eachindex(shifts), station in 1:3
        a, b = rand(rng, 50), rand(rng, 50)
        table[k, station, 1, :] = [50.0, sum(a), sum(b), sum(abs2, a), sum(abs2, b), sum(a .* b)]
    end
    table[2, :, 1, :] = table[8, :, 1, :]    # (0,-1) and (0,+1) tie exactly
    train = [1, 3]
    pooled = [fm_pooled_r(vec(sum(table[k, train, 1, :]; dims=1))) for k in eachindex(shifts)]
    order = sortperm(shifts; by=s -> (s[1]^2 + s[2]^2, s[2], s[1]))
    best = order[1]
    for k in order
        pooled[k] > pooled[best] && (best = k)
    end
    k, r, r0 = FM_SC.choose_shift(grid, table, train, 1)
    @test k == best
    @test r ≈ pooled[best]
    @test r0 ≈ pooled[findfirst(==((0, 0)), shifts)]
end

"""The wet threshold and the 199 knots of one cell, from its training pairs."""
function fm_fit_map(x, y; wet=0.1)
    w = mean(y .>= wet)
    tau = quantile(x, 1 - w)
    probabilities = collect(0.005:0.005:0.995)
    return (; tau, qx=quantile(x[x .> tau], probabilities), qy=quantile(y[y .>= wet], probabilities))
end

function fm_map(m, x)
    x <= m.tau && return 0.0
    x >= m.qx[end] && return m.qy[end] * x / m.qx[end]
    x <= m.qx[1] && return m.qy[1]
    j = findlast(<=(x), m.qx)
    return m.qy[j] + (x - m.qx[j]) / (m.qx[j + 1] - m.qx[j]) * (m.qy[j + 1] - m.qy[j])
end

@testset "formula: wet threshold and quantile map" begin
    rng = MersenneTwister(3)
    n_station, n_time = 12, 900
    rain(scale) = [rand(rng) < 0.25 ? scale * randexp(rng) : 0.0 for _ in 1:n_station, _ in 1:n_time]
    y_obs, y_sat = round.(rain(2.0) .* 2) ./ 2, rain(3.0) .+ 0.02 .* rand(rng, n_station, n_time)
    train = collect(1:9)
    cells = ones(Int, n_time)
    maps = FM_SC.fit_quantile_map(y_obs, y_sat, train, cells; min_wet=50)
    ours = fm_fit_map(vec(y_sat[train, :]), vec(y_obs[train, :]))
    @test maps[1].tau ≈ ours.tau
    @test maps[1].qx ≈ ours.qx
    @test maps[1].qy ≈ ours.qy
    mapped = FM_SC.apply_quantile_map(maps, y_sat, cells)
    @test fm_close(mapped, [fm_map(ours, y_sat[i, t]) for i in 1:n_station, t in 1:n_time])
    @test fm_map(ours, 2 * ours.qx[end]) ≈ 2 * ours.qy[end]         # the upper tail scales
    @test count(mapped[train, :] .> 0) == count(y_sat[train, :] .> ours.tau)

    # A sparse cell takes its regime's map: cells 1 and 2 are seasons of regime 1.
    cells2 = [t <= 880 ? 1 : 2 for t in 1:n_time]
    maps2 = FM_SC.fit_quantile_map(y_obs, y_sat, train, cells2; min_wet=200)
    @test maps2[1].pool == "cell"
    @test maps2[2].pool == "regime"
    @test maps2[2].tau ≈ ours.tau
end

"""Candidate c of the adaptive strength: 0 leaves x alone, 1..5 are the threshold with 0, ¼, ½, ¾, 1 of the map."""
function fm_adjust(x, mapped, tau, c)
    c == 0 && return x
    x <= tau && return 0.0
    s = (0.0, 0.25, 0.5, 0.75, 1.0)[c]
    return (1 - s) * x + s * mapped
end

@testset "formula: strength of the map" begin
    rng = MersenneTwister(4)
    n_station, n_time = 10, 600
    y_obs = [rand(rng) < 0.3 ? randexp(rng) : 0.0 for _ in 1:n_station, _ in 1:n_time]
    y_sat = [rand(rng) < 0.3 ? 2 * randexp(rng) : 0.0 for _ in 1:n_station, _ in 1:n_time]
    y_sat[:, 1:200] .= y_obs[:, 1:200] .* 1.5       # month 1: the map helps
    train = collect(1:8)
    cells = ones(Int, n_time)
    months = [t <= 200 ? 1 : t <= 400 ? 2 : 3 for t in 1:n_time]
    maps = FM_SC.fit_quantile_map(y_obs, y_sat, train, cells; min_wet=50)
    m = fm_fit_map(vec(y_sat[train, :]), vec(y_obs[train, :]))
    chosen = FM_SC.choose_strength(y_obs, y_sat, maps, train, cells, months)
    for month in 1:3
        sse = [sum((fm_adjust(y_sat[i, t], fm_map(m, y_sat[i, t]), m.tau, c) - y_obs[i, t])^2
                   for i in train, t in findall(==(month), months)) for c in 0:5]
        @test chosen[month].candidate == argmin(sse) - 1
        @test chosen[month].rmse_none ≈ sqrt(sse[1] / (length(train) * count(==(month), months)))
    end
    @test FM_SC.QM_STRENGTHS == (0.0, 0.25, 0.5, 0.75, 1.0)
end

# --- C / E. Distances, interpolators -------------------------------------------------------------

function fm_haversine(lon1, lat1, lon2, lat2)
    p1, p2 = lat1 * pi / 180, lat2 * pi / 180
    a = sin((p2 - p1) / 2)^2 + cos(p1) * cos(p2) * sin((lon2 - lon1) * pi / 360)^2
    return 2 * FM_R_KM * asin(sqrt(a))
end

fm_local_km(lon, lat, lon0, lat0) =
    ((lon - lon0) * cosd(lat0) * pi / 180 * FM_R_KM, (lat - lat0) * pi / 180 * FM_R_KM)

function fm_stations(rng, n)
    return hcat(109.5 .+ 2 .* rand(rng, n), 31.5 .+ 1.5 .* rand(rng, n))
end

@testset "formula: distances and local coordinates" begin
    rng = MersenneTwister(5)
    a, b = fm_stations(rng, 6), fm_stations(rng, 4)
    D = haversine_distance_matrix(a, b)
    @test fm_close(D, [fm_haversine(a[i, 1], a[i, 2], b[j, 1], b[j, 2]) for i in 1:6, j in 1:4])
    xy = FM_TI.local_km_coordinates(a)
    lon0, lat0 = mean(a[:, 1]), mean(a[:, 2])
    @test fm_close(xy, reduce(vcat, [collect(fm_local_km(a[i, 1], a[i, 2], lon0, lat0))' for i in 1:6]))
end

"""IDW or ADW at one target and hour, over the reporting stations (the `k` nearest of them, or all)."""
function fm_idw(train, values, target; power, k=nothing, angular=false)
    reporting = findall(isfinite, values)
    d = [fm_haversine(train[i, 1], train[i, 2], target[1], target[2]) for i in reporting]
    if k !== nothing
        keep = sortperm(d)[1:min(k, length(d))]
        reporting, d = reporting[keep], d[keep]
    end
    on_top = findall(<=(sqrt(eps())), d)
    isempty(on_top) || return mean(values[reporting[on_top]])
    w = d .^ (-power)
    if angular && length(reporting) > 1
        u = [collect(fm_local_km(train[i, 1], train[i, 2], target[1], target[2])) for i in reporting]
        u = [v ./ hypot(v...) for v in u]
        w = [w[j] * (1 + sum(w[l] * (1 - dot(u[j], u[l])) for l in eachindex(w) if l != j) /
                         sum(w[l] for l in eachindex(w) if l != j)) for j in eachindex(w)]
    end
    return sum(w .* values[reporting]) / sum(w)
end

@testset "formula: IDW and ADW" begin
    rng = MersenneTwister(6)
    train, target = fm_stations(rng, 14), fm_stations(rng, 5)
    values = rand(rng, 14, 6) .* 5
    values[3, 2] = NaN; values[7, 2] = NaN; values[1, 5] = NaN
    for power in (1.0, 2.0, 3.0), k in (nothing, 4)
        idw = FM_TI.idw_predict(train, values, target; power, neighbors=k)
        adw = FM_TI.adw_predict(train, values, target; power, neighbors=k)
        @test fm_close(idw, [fm_idw(train, values[:, t], target[j, :]; power, k) for j in 1:5, t in 1:6]; tol=1e-9)
        @test fm_close(adw, [fm_idw(train, values[:, t], target[j, :]; power, k, angular=true) for j in 1:5, t in 1:6]; tol=1e-9)
    end
    # A gauge on top of the target is returned as it is.
    @test FM_TI.adw_predict(train, values, train[4:4, :]; power=2.0)[1, 1] ≈ values[4, 1]
end

fm_tps_kernel(r2) = r2 <= eps() ? 0.0 : 0.5 * r2 * log(r2)

"""Thin-plate smoothing spline at the targets for one hour; `scale` is fixed over the whole call."""
function fm_tps(xy, y, target_xy, smooth, scale)
    n = size(xy, 1)
    K = [fm_tps_kernel(sum(abs2, xy[i, :] .- xy[j, :])) for i in 1:n, j in 1:n]
    P = hcat(ones(n), xy)
    coefficients = [K + smooth * scale * I P; P' zeros(3, 3)] \ vcat(y, zeros(3))
    a, b = coefficients[1:n], coefficients[(n + 1):end]
    return [sum(a[j] * fm_tps_kernel(sum(abs2, target_xy[t, :] .- xy[j, :])) for j in 1:n) +
            b[1] + b[2] * target_xy[t, 1] + b[3] * target_xy[t, 2] for t in axes(target_xy, 1)]
end

@testset "formula: thin-plate spline" begin
    rng = MersenneTwister(7)
    train, target = fm_stations(rng, 12), fm_stations(rng, 4)
    values = rand(rng, 12, 3) .* 4
    values[5, 2] = NaN
    lon0, lat0 = mean(train[:, 1]), mean(train[:, 2])
    xy = reduce(vcat, [collect(fm_local_km(train[i, 1], train[i, 2], lon0, lat0))' for i in 1:12])
    txy = reduce(vcat, [collect(fm_local_km(target[i, 1], target[i, 2], lon0, lat0))' for i in 1:4])
    K = [fm_tps_kernel(sum(abs2, xy[i, :] .- xy[j, :])) for i in 1:12, j in 1:12]
    scale = median(abs.(K[K .!= 0]))
    for smooth in (1e-3, 0.1, 1.0)
        theirs = FM_TI.tps_predict(train, values, target; smooth)
        for t in 1:3
            keep = findall(isfinite, values[:, t])
            @test fm_close(theirs[:, t], fm_tps(xy[keep, :], values[keep, t], txy, smooth, scale); tol=1e-7)
        end
    end
end

# --- E. Kernels, GWR, mixed GWR, MGWR ------------------------------------------------------------

fm_kernel(name, d, b) = name == GAUSSIAN ? exp(-d^2 / (2 * b^2)) :
    name == EXPONENTIAL ? exp(-d / b) :
    name == BISQUARE ? (d > b ? 0.0 : (1 - (d / b)^2)^2) :
    name == TRICUBE ? (d > b ? 0.0 : (1 - (d / b)^3)^3) : (d > b ? 0.0 : 1.0)

"""Adaptive bandwidth: the distance to the k-th nearest station, or k/m of the largest when only m < k exist."""
fm_bandwidth(d, k) = k <= length(d) ? sort(d)[k] : k / length(d) * maximum(d)

@testset "formula: kernels and the adaptive bandwidth" begin
    rng = MersenneTwister(8)
    d = 100 .* rand(rng, 20)
    for kernel in (GAUSSIAN, EXPONENTIAL, BISQUARE, TRICUBE, BOXCAR)
        @test gw_weight(d, 35.0; kernel, adaptive=false) ≈ [fm_kernel(kernel, x, 35.0) for x in d]
        for k in (5, 20, 30)
            @test gw_weight(d, Float64(k); kernel, adaptive=true) ≈ [fm_kernel(kernel, x, fm_bandwidth(d, k)) for x in d]
        end
        @test all(gw_weight(d, Inf; kernel, adaptive=false) .== 1.0)   # an infinite bandwidth is a global fit
    end
end

"""Row of the local hat matrix at one target: x0' (X'WX + ridge I)^-1 X'W."""
fm_hat_row(X, x0, w, ridge) = vec(x0' * ((X' * (w .* X) + ridge * I) \ (X' .* w')))

function fm_hat(X, Xt, D, k, kernel, ridge; exclude_self=false)
    H = zeros(size(Xt, 1), size(X, 1))
    for t in axes(Xt, 1)
        d = D[:, t]
        others = exclude_self ? setdiff(1:length(d), t) : collect(1:length(d))
        b = fm_bandwidth(d[others], k)
        w = [i in others ? fm_kernel(kernel, d[i], b) : 0.0 for i in eachindex(d)]
        H[t, :] = fm_hat_row(X, Xt[t, :], w, ridge)
    end
    return H
end

@testset "formula: local fit (hat matrix)" begin
    rng = MersenneTwister(9)
    train, target = fm_stations(rng, 25), fm_stations(rng, 6)
    X, Xt = hcat(ones(25), randn(rng, 25, 2)), hcat(ones(6), randn(rng, 6, 2))
    D = haversine_distance_matrix(train, target)
    for kernel in (GAUSSIAN, EXPONENTIAL, BISQUARE)
        H = local_hat_matrix(X, Xt, D, 12.0, kernel; adaptive=true, ridge=1e-8)
        @test fm_close(H, fm_hat(X, Xt, D, 12, kernel, 1e-8); tol=1e-8)
    end
    # Fewer than p + 1 stations with weight: the target is not fitted.
    H = local_hat_matrix(X, Xt, D, 1.0, BISQUARE; adaptive=false, unsupported=:missing)
    @test all(isnan, H)
end

@testset "formula: mixed GWR back-fit" begin
    rng = MersenneTwister(10)
    train, target = fm_stations(rng, 30), fm_stations(rng, 5)
    Xl, Xlt = hcat(ones(30), randn(rng, 30)), hcat(ones(5), randn(rng, 5))
    Xg, Xgt = randn(rng, 30, 1), randn(rng, 5, 1)
    y = randn(rng, 30, 1)
    ridge, tol = 1e-8, 1e-5
    theirs, ok = mixed_gwr_predict(Xl, Xg, y, train, Xlt, Xgt, target, 12.0, EXPONENTIAL;
        adaptive=true, ridge, tolerance=tol, max_iterations=1000)
    @test all(ok)
    H = fm_hat(Xl, Xl, haversine_distance_matrix(train, train), 12, EXPONENTIAL, ridge)
    P = Xg * ((Xg' * Xg + ridge * I) \ Xg')
    v = vec(y)
    L, g, previous = zeros(30), P * v, Inf
    for _ in 1:1000
        L = H * (v - g)
        g = P * (v - L)
        rss = sum(abs2, v - L - g)
        change = isfinite(previous) ? abs(previous - rss) / max(previous, eps()) : Inf
        previous = rss
        change < tol && break
    end
    beta = (Xg' * Xg + ridge * I) \ (Xg' * (v - L))
    Ht = fm_hat(Xl, Xlt, haversine_distance_matrix(train, target), 12, EXPONENTIAL, ridge)
    @test fm_close(vec(theirs), Ht * (v - Xg * beta) + Xgt * beta; tol=1e-8)
end

@testset "formula: MGWR back-fit with over-relaxation" begin
    rng = MersenneTwister(11)
    train, target = fm_stations(rng, 40), fm_stations(rng, 6)
    groups = [ones(40, 1), randn(rng, 40, 1), randn(rng, 40, 1)]
    groups_t = [ones(6, 1), randn(rng, 6, 1), randn(rng, 6, 1)]
    bandwidths = [8.0, 20.0, 30.0]
    y = randn(rng, 40, 1)
    cfg = (; ridge=1e-8, relaxation=1.5, max_iterations=1000, tolerance=1e-5)
    for (Xg, Xgt) in ((zeros(40, 0), zeros(6, 0)), (randn(rng, 40, 1), randn(rng, 6, 1)))
        theirs, ok = FM_JCM._multiscale_predict_damped(groups, Xg, y, train, groups_t, Xgt, target,
            bandwidths, EXPONENTIAL, cfg; adaptive=true)
        @test all(ok)
        Dtt, Dt = haversine_distance_matrix(train, train), haversine_distance_matrix(train, target)
        H = [fm_hat(groups[j], groups[j], Dtt, Int(bandwidths[j]), EXPONENTIAL, cfg.ridge) for j in 1:3]
        has_global = !isempty(Xg)
        P = has_global ? Xg * ((Xg' * Xg + cfg.ridge * I) \ Xg') : zeros(40, 40)
        v, omega = vec(y), cfg.relaxation
        f = [zeros(40) for _ in 1:3]
        g = has_global ? P * v : zeros(40)
        previous = fill(NaN, 40)
        for _ in 1:cfg.max_iterations
            for j in 1:3
                f[j] = omega .* (H[j] * (v - g - sum(f[l] for l in 1:3 if l != j))) .+ (1 - omega) .* f[j]
            end
            has_global && (g = omega .* (P * (v - sum(f))) .+ (1 - omega) .* g)
            fitted = sum(f) + g
            change = all(isfinite, previous) ? norm(fitted - previous) / max(norm(previous), eps()) : Inf
            previous = fitted
            change < cfg.tolerance && break
        end
        beta = has_global ? (Xg' * Xg + cfg.ridge * I) \ (Xg' * (v - sum(f))) : Float64[]
        prediction = has_global ? Xgt * beta : zeros(6)
        for j in 1:3
            partial = v - (has_global ? Xg * beta : zeros(40)) - sum(f[l] for l in 1:3 if l != j)
            prediction += fm_hat(groups[j], groups_t[j], Dt, Int(bandwidths[j]), EXPONENTIAL, cfg.ridge) * partial
        end
        @test fm_close(vec(theirs), prediction; tol=1e-7)
    end
end

@testset "formula: covariate centring and scaling" begin
    rng = MersenneTwister(12)
    n, nt = 12, 40
    lonlat = fm_stations(rng, n)
    y_obs = [rand(rng) < 0.4 ? randexp(rng) : 0.0 for _ in 1:n, _ in 1:nt]
    y_sat = rand(rng, n, nt)
    covariate = 10 .+ randn(rng, n, nt) .+ (1:nt)' ./ 5
    train, target = collect(1:9), [10, 11, 12]
    cfg = FM_JCM.JointCovariateBenchmarkConfig(terrain_path="", era5_paths=Dict{Int,String}())
    context = FM_JCM.build_joint_fold_context("TEST", Dict("t2m_c" => "local"), train, target, lonlat,
        y_obs, y_sat, DataFrame(), Dict(:t2m_c => covariate), nothing, cfg)
    # Each hour loses its training mean; then one centre and scale over the training wet hours,
    # every station weighing the same in total.
    centred = covariate .- mean(covariate[train, :]; dims=1)
    wet = y_obs[train, :] .>= 0.1
    weights = wet ./ max.(sum(wet; dims=2), 1)
    centre = sum(weights .* centred[train, :]) / sum(weights)
    scale = sqrt(sum(weights .* (centred[train, :] .- centre) .^ 2) / sum(weights))
    @test fm_close(context.predictor_train["t2m_c"][:, 1, :], (centred[train, :] .- centre) ./ scale; tol=1e-9)
    @test fm_close(context.predictor_target["t2m_c"][:, 1, :], (centred[target, :] .- centre) ./ scale; tol=1e-9)
    # The spatial block: 1 and the coordinates standardised over the training gauges.
    z(col) = (lonlat[:, col] .- mean(lonlat[train, col])) ./ std(lonlat[train, col])
    @test fm_close(context.spatial_train, hcat(ones(9), z(1)[train], z(2)[train]))
end

# --- E. Covariate selection statistics ------------------------------------------------------------

@testset "formula: VIF, permutation p, Benjamini-Hochberg, joint F" begin
    rng = MersenneTwister(13)
    X = randn(rng, 200, 3)
    X[:, 3] .+= 0.8 .* X[:, 1]
    w = rand(rng, 200)
    function vif(j)
        others = setdiff(1:3, j)
        design = hcat(ones(200), X[:, others])
        beta = (design' * (w .* design)) \ (design' * (w .* X[:, j]))
        m = sum(w .* X[:, j]) / sum(w)
        r2 = 1 - sum(w .* (X[:, j] .- design * beta) .^ 2) / sum(w .* (X[:, j] .- m) .^ 2)
        return 1 / (1 - r2)
    end
    @test FM_JVS._column_vifs(X, w) ≈ [vif(j) for j in 1:3]

    p = rand(rng, 9)
    order = sortperm(p)
    q = similar(p)
    for (rank, index) in enumerate(order)
        q[index] = minimum(min(1.0, p[order[j]] * 9 / j) for j in rank:9)
    end
    @test FM_DEM.bh_adjust(p) ≈ q

    x, y = randn(rng, 60, 1), randn(rng, 60)
    statistic(a, b) = abs(cor(vec(a), b))
    observed = statistic(x, y)
    theirs = FM_DEM._permutation_pvalue(statistic, x, y, observed, 199, MersenneTwister(99))
    shuffled, exceed, rng2 = copy(y), 0, MersenneTwister(99)
    for _ in 1:199
        shuffle!(rng2, shuffled)
        exceed += statistic(x, shuffled) >= observed
    end
    @test theirs == (exceed + 1) / 200

    A = randn(rng, 60, 2)
    design = hcat(ones(60), A)
    rss1 = sum(abs2, y - design * (design \ y))
    rss0 = sum(abs2, y .- mean(y))
    @test FM_DEM._joint_f_statistic(A, y) ≈ ((rss0 - rss1) / 2) / (rss1 / (60 - 2 - 1)) rtol = 1e-6
end

# --- E. Tuning hours, merges, blends ----------------------------------------------------------------

@testset "formula: tuning hours" begin
    rng = MersenneTwister(14)
    y = rand(rng, 8, 1000)
    hours = _tuning_time_indices(y, 336)
    wettest = sortperm(vec(mean(y; dims=1)); rev=true)[1:168]
    spaced = unique(round.(Int, range(1, 1000; length=168)))
    @test hours == sort(unique(vcat(wettest, spaced)))
    @test _tuning_time_indices(y[:, 1:300], 336) == collect(1:300)
end

@testset "formula: merged anchors" begin
    rng = MersenneTwister(15)
    sources = [rand(rng, 10, 50) .* 3 for _ in 1:3]
    y_obs = 0.1 .+ 0.2 .* sources[1] .+ 0.5 .* sources[2] .+ 0.1 .* sources[3] .+ 0.3 .* randn(rng, 10, 50)
    y_obs[2, 7] = NaN; sources[2][4, 9] = NaN
    train = collect(1:7)
    beta, used, fell_back = FM_SF.fusion_coefficients(:ols, y_obs, sources, train)
    keep = [(i, t) for i in train, t in 1:50 if isfinite(y_obs[i, t]) && all(isfinite(s[i, t]) for s in sources)]
    design = [k == 0 ? 1.0 : sources[k][i, t] for (i, t) in keep, k in 0:3]
    @test !fell_back && used == length(keep)
    @test beta ≈ design \ [y_obs[i, t] for (i, t) in keep]
    merged = FM_SF.apply_satellite_fusion(sources, beta)
    @test isnan(merged[4, 9])
    @test merged[1, 1] ≈ max(beta[1] + sum(beta[k + 1] * sources[k][1, 1] for k in 1:3), 0.0)
    equal, _, _ = FM_SF.fusion_coefficients(:mean, y_obs, sources, train)
    @test FM_SF.apply_satellite_fusion(sources, equal)[3, 3] ≈ mean(s[3, 3] for s in sources)
end

"""Agreement band of one cell: 0 when no product is wet, else 3(a - 1) + e."""
function fm_band(values; wet=0.1, edges=(0.1, 2.5, 8.0))
    finite = filter(isfinite, values)
    a = min(count(>=(wet), finite), 3)
    a == 0 && return 0
    return 3 * (a - 1) + count(<=(maximum(finite)), edges)
end

@testset "formula: blends and agreement bands" begin
    rng = MersenneTwister(16)
    shape = (9, 400)
    field() = [rand(rng) < 0.35 ? 6 * randexp(rng) : 0.0 for _ in 1:shape[1], _ in 1:shape[2]]
    sources = [field() for _ in 1:3]
    sources[2][1, 1] = NaN
    y_obs, prediction, fallback = field(), field(), field()
    y_sat = sources[1]
    band = blend_band_matrix(:agreement_envelope, y_sat, sources, 0.1)
    @test band == [fm_band([s[i, t] for s in sources]) for i in 1:shape[1], t in 1:shape[2]]
    @test extrema(band) == (0, 9)

    blended = satellite_wet_blend_prediction(y_sat, prediction, fallback, 0.3, 0.1)
    @test fm_close(blended, [y_sat[i] >= 0.1 ? 0.7 * prediction[i] + 0.3 * fallback[i] : prediction[i]
                             for i in CartesianIndices(prediction)])

    lambdas = collect(0.0:0.1:1.0)
    selected = select_blend_lambdas(y_obs, y_sat, band, prediction, fallback, lambdas, 9)
    for b in 1:9
        cells = findall(==(b), band)
        sse = [sum(((1 - l) * prediction[c] + l * fallback[c] - y_obs[c])^2 for c in cells) for l in lambdas]
        @test selected.lambdas[b] == lambdas[argmin(sse)]
    end
    banded = banded_blend_prediction(band, prediction, fallback, selected.lambdas)
    @test banded[band .== 0] == prediction[band .== 0]
end

# --- F. Scores and the bootstrap -----------------------------------------------------------------------

@testset "formula: scores" begin
    rng = MersenneTwister(17)
    obs = [rand(rng) < 0.3 ? 3 * randexp(rng) : 0.0 for _ in 1:2000]
    est = [rand(rng) < 0.3 ? 3 * randexp(rng) : 0.0 for _ in 1:2000] .+ 0.4 .* obs
    c = metric_continuous(obs, est)
    @test c.RMSE ≈ sqrt(mean((est .- obs) .^ 2))
    @test c.MAE ≈ mean(abs.(est .- obs))
    @test c.Bias ≈ mean(est .- obs)
    @test c.r ≈ cor(obs, est)
    for threshold in (0.1, 2.0, 8.0)
        hits = count((obs .>= threshold) .& (est .>= threshold))
        misses = count((obs .>= threshold) .& (est .< threshold))
        false_alarms = count((obs .< threshold) .& (est .>= threshold))
        e = metric_event(obs, est; thr=threshold)
        @test e.POD ≈ hits / (hits + misses)
        @test e.FAR ≈ false_alarms / (hits + false_alarms)
        @test e.CSI ≈ hits / (hits + misses + false_alarms)
    end
end

@testset "formula: paired day-block bootstrap and Holm" begin
    rng = MersenneTwister(18)
    times = collect(DateTime(2022, 6, 1):Hour(1):DateTime(2022, 6, 20, 23))
    shape = (6, length(times))
    y_obs, baseline, method = rand(rng, shape...), rand(rng, shape...), rand(rng, shape...)
    mask = BitMatrix(rand(rng, shape...) .< 0.9)
    mask[:, 25:48] .= false                              # a day with no cell is never drawn
    theirs = _daily_bootstrap_delta(MersenneTwister(7), times, y_obs, baseline, method, mask, 300)
    days = Date.(times)
    blocks = [findall(==(day), days) for day in sort(unique(days))]
    n = [count(mask[:, b]) for b in blocks]
    sse(prediction) = [sum(abs2, (prediction[:, b] .- y_obs[:, b])[mask[:, b]]) for b in blocks]
    sse_base, sse_method = sse(baseline), sse(method)
    keep = findall(>(0), n)
    rng2 = MersenneTwister(7)
    ours = map(1:300) do _
        drawn = rand(rng2, keep, length(keep))
        sqrt(sum(sse_base[drawn]) / sum(n[drawn])) - sqrt(sum(sse_method[drawn]) / sum(n[drawn]))
    end
    @test theirs ≈ ours

    p = [0.01, 0.04, 0.03, 0.2]
    sorted = sort(p)
    adjusted = [maximum(min(1.0, (4 - j + 1) * sorted[j]) for j in 1:findfirst(==(value), sorted)) for value in p]
    @test _holm_adjust(p) ≈ adjusted
end
