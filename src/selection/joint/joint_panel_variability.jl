function _panel_stats(panel, groups::Vector{String}, lonlat::Matrix{Float64}, cfg)
    eligible = findall(panel.eligible_station)
    isempty(eligible) && return nothing
    coords = lonlat[eligible, :]
    means = vec(mean(coords, dims=1)); scales = vec(std(coords, dims=1))
    any(x -> !isfinite(x) || x <= 0, scales) && return nothing
    columns = reduce(vcat, [GROUP_COLUMNS[group] for group in groups]; init=Symbol[])
    p = 3 + length(columns)
    A = [zeros(p, p) for _ in eligible]
    b = [zeros(p) for _ in eligible]
    rows = [Vector{Vector{Float64}}() for _ in eligible]
    ys = [Float64[] for _ in eligible]
    row_weights = [Float64[] for _ in eligible]
    for (position, i) in enumerate(eligible), time in axes(panel.y, 2)
        weight = panel.weights[i, time]
        weight > 0 && isfinite(panel.y[i, time]) &&
            all(isfinite(panel.x[column][i, time]) for column in columns) || continue
        design = [1.0, (coords[position, 1] - means[1]) / scales[1],
            (coords[position, 2] - means[2]) / scales[2],
            (panel.x[column][i, time] for column in columns)...]
        A[position] .+= weight .* (design * design')
        b[position] .+= weight .* design .* panel.y[i, time]
        push!(rows[position], design); push!(ys[position], panel.y[i, time])
        push!(row_weights[position], weight)
    end
    group_indices = Dict{String,Vector{Int}}()
    offset = 3
    for group in groups
        width = length(GROUP_COLUMNS[group])
        group_indices[group] = collect((offset + 1):(offset + width))
        offset += width
    end
    return (; A, b, rows, ys, row_weights, coords, p, group_indices, eligible)
end

function _coefficients(stats, distances, neighbors, cfg;
    order=collect(eachindex(stats.A)), weight_rows=nothing)
    n = length(stats.A); result = fill(NaN, n, stats.p)
    # Weights are a function of the geometry and bandwidth alone - `order` reorders which
    # sufficient statistic each spatial slot contributes, never the weight on that slot - so a
    # permutation test builds them once and hands them in.
    rows = weight_rows === nothing ?
        eachcol(gw_weight(distances, Float64(neighbors); kernel=BISQUARE)) : weight_rows
    lhs = Matrix{Float64}(undef, stats.p, stats.p); rhs = Vector{Float64}(undef, stats.p)
    for target in 1:n
        spatial_weights = rows[target]
        # Refilled to exactly what `cfg.ridge .* Matrix{Float64}(I, p, p)` / `zeros(p)` produced,
        # rather than reallocated per target.
        fill!(lhs, 0.0); fill!(rhs, 0.0)
        for i in 1:stats.p
            lhs[i, i] = cfg.ridge
        end
        for spatial_index in 1:n
            source = order[spatial_index]
            lhs .+= spatial_weights[spatial_index] .* stats.A[source]
            rhs .+= spatial_weights[spatial_index] .* stats.b[source]
        end
        try result[target, :] .= lhs \ rhs catch end
    end
    return result
end

function _loocv(stats, distances, neighbors, cfg)
    total = 0.0; station_count = 0; n = length(stats.A)
    loo = loocv_distances(distances)
    for target in 1:n
        spatial_weights = gw_weight(view(loo, :, target), Float64(neighbors); kernel=BISQUARE)
        lhs = cfg.ridge .* Matrix{Float64}(I, stats.p, stats.p); rhs = zeros(stats.p)
        for i in 1:n
            lhs .+= spatial_weights[i] .* stats.A[i]; rhs .+= spatial_weights[i] .* stats.b[i]
        end
        beta = try lhs \ rhs catch; continue end
        isempty(stats.rows[target]) && continue
        station_sse = sum(stats.row_weights[target][j] *
            (stats.ys[target][j] - dot(stats.rows[target][j], beta))^2
            for j in eachindex(stats.ys[target])) / sum(stats.row_weights[target])
        total += station_sse; station_count += 1
    end
    return station_count == n ? sqrt(total / n) : Inf
end

function joint_spatial_variability_test(
    panel, groups::Vector{String}, lonlat::Matrix{Float64}, cfg::JointSelectionConfig;
    rng::AbstractRNG=MersenneTwister(cfg.seed),
)
    scan = DataFrame(neighbors=Int[], loocv_rmse=Float64[], available=Bool[])
    result = DataFrame(
        variable_group=String[], statistic=Float64[], pvalue=Float64[], qvalue=Float64[],
        role=String[], neighbors=Union{Missing,Int}[], status=String[],
    )
    isempty(groups) && return (; scan, result, bandwidth=missing)
    stats = _panel_stats(panel, groups, lonlat, cfg)
    if stats === nothing
        for group in groups
            push!(result, (group, NaN, NaN, NaN, "uncertain", missing, "insufficient_stations"))
        end
        return (; scan, result, bandwidth=missing)
    end
    distances = haversine_distance_matrix(stats.coords, stats.coords)
    candidates = filter(c -> max(4, stats.p + 1) <= c <= length(stats.eligible) - 1,
        unique(cfg.bandwidth_candidates))
    for candidate in candidates
        value = _loocv(stats, distances, candidate, cfg)
        push!(scan, (candidate, value, isfinite(value)))
    end
    available = filter(:available => identity, scan)
    if nrow(available) == 0
        for group in groups
            push!(result, (group, NaN, NaN, NaN, "uncertain", missing, "no_bandwidth"))
        end
        return (; scan, result, bandwidth=missing)
    end
    bandwidth = available.neighbors[argmin(available.loocv_rmse)]
    weight_rows = eachcol(gw_weight(distances, Float64(bandwidth); kernel=BISQUARE))
    coefficients = _coefficients(stats, distances, bandwidth, cfg; weight_rows)
    observed = [sum(var(coefficients[:, index]) for index in stats.group_indices[group]) for group in groups]
    # Drawn up front, from the same `rng` in the same order, so `orders[i]` is exactly what
    # iteration `i` drew for itself; the bodies are then independent and each is a full
    # O(n^2 p^2) coefficient-surface rebuild. `hits` is one row per permutation, reduced
    # afterwards - summing booleans is exact and order-free, so `exceed` is unchanged.
    orders = [station_block_permutation(length(stats.A), rng)
        for _ in 1:cfg.spatial_permutations]
    # `Matrix{Bool}`, not a `BitMatrix`: a BitArray packs 64 entries into one word, so two threads
    # writing different rows of the same column would read-modify-write the same word and lose
    # updates. One byte per entry makes the concurrent writes independent.
    hits = fill(false, cfg.spatial_permutations, length(groups))
    Threads.@threads :greedy for permutation_index in 1:cfg.spatial_permutations
        permuted = _coefficients(stats, distances, bandwidth, cfg;
            order=orders[permutation_index], weight_rows)
        for (group_index, group) in enumerate(groups)
            value = sum(var(permuted[:, index]) for index in stats.group_indices[group])
            hits[permutation_index, group_index] = isfinite(value) && value >= observed[group_index]
        end
    end
    exceed = vec(sum(hits, dims=1))
    pvalues = (exceed .+ 1) ./ (cfg.spatial_permutations + 1)
    qvalues = DEMTerrainExperiment.bh_adjust(pvalues)
    for (index, group) in enumerate(groups)
        role = qvalues[index] < cfg.q_threshold ? "local" : "global"
        push!(result, (group, observed[index], pvalues[index], qvalues[index],
            role, bandwidth, "ok"))
    end
    return (; scan, result, bandwidth)
end
