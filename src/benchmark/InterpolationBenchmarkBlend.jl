"""
Blend a residual-family prediction toward a gauge-only one on the cells the satellite calls wet.

`lambda = 0` returns `prediction` untouched, `lambda = 1` replaces it with `fallback` on those
cells and leaves the satellite-dry cells alone either way. The test is on the satellite, which is
knowable at prediction time - this is a model, not an oracle.

A NaN in `prediction` or the satellite propagates. A NaN in `fallback` does not: the cell simply
keeps its unblended prediction, so a fallback that failed on some cells cannot shrink the blended
method's coverage below its source's. That matters because coverage feeds the run's status rows.
"""
function satellite_wet_blend_prediction(
    y_sat::Matrix{Float64}, prediction::Matrix{Float64}, fallback::Matrix{Float64},
    lambda::Float64, threshold::Float64,
)
    size(prediction) == size(y_sat) == size(fallback) || throw(DimensionMismatch(
        "satellite, prediction and fallback must have the same shape",
    ))
    out = similar(prediction)
    @inbounds for index in eachindex(prediction)
        satellite = y_sat[index]
        predicted = prediction[index]
        other = fallback[index]
        if isnan(satellite) || isnan(predicted)
            out[index] = NaN
        elseif satellite >= threshold && !isnan(other)
            out[index] = (1 - lambda) * predicted + lambda * other
        else
            out[index] = predicted
        end
    end
    return out
end

"""
Blend toward the fallback with one weight per band, the band read off a precomputed matrix.

The vector-weight twin of `satellite_wet_blend_prediction`, and it keeps that function's rules
exactly: band 0 is never touched, a NaN prediction propagates, and a NaN fallback leaves the cell
at its unblended value rather than poisoning it. That last one is what makes the per-band
selection in `select_blend_lambdas` exact - the set of finite cells does not move as the weights
move, so pooled SSE is additive over the bands and each band minimises independently.

`band` is `0` for a cell no weight reaches, and otherwise an index into `lambdas`.
"""
function banded_blend_prediction(
    band::Matrix{Int}, prediction::Matrix{Float64}, fallback::Matrix{Float64},
    lambdas::Vector{Float64},
)
    size(prediction) == size(band) == size(fallback) || throw(DimensionMismatch(
        "band, prediction and fallback must have the same shape",
    ))
    out = similar(prediction)
    @inbounds for index in eachindex(prediction)
        predicted = prediction[index]
        other = fallback[index]
        position = band[index]
        if isnan(predicted)
            out[index] = NaN
        elseif position == 0 || isnan(other)
            out[index] = predicted
        else
            1 <= position <= length(lambdas) || throw(BoundsError(lambdas, position))
            lambda = lambdas[position]
            out[index] = (1 - lambda) * predicted + lambda * other
        end
    end
    return out
end

"""
Band index per cell for `axis`, from the source products the fusion also draws on.

`:constant` reproduces the scalar intervention set exactly - band 1 wherever the *target's* own
anchor calls the cell wet - so a one-band run of the banded path is the scalar path.

`:agreement_envelope` bands on the source products instead: `(agreement - 1) * 3 + envelope`,
where `agreement` counts how many products call the cell wet and `envelope` places `max` over them
in `BLEND_ENVELOPE_EDGES`. A NaN source counts as not wet rather than poisoning the cell, so a
cell can still be reached when one product is missing and the count is then a lower bound - which
biases the weight toward the unblended model rather than toward blending on absent evidence. Any
product being wet puts the envelope at or above the first edge, so a nonzero agreement always
pairs with a nonzero envelope band.
"""
function blend_band_matrix(
    axis::Symbol, y_sat::Matrix{Float64}, sources::Vector{Matrix{Float64}}, threshold::Float64,
)
    axis in BLEND_AXES || throw(ArgumentError("unknown blend axis: $axis"))
    if axis in BLEND_DRY_AXES
        # The agreement axis's bands, with its never-blended band 0 given the last ID instead.
        out = blend_band_matrix(:agreement_envelope, y_sat, sources, threshold)
        out[out .== 0] .= blend_band_count(axis)
        return out
    end
    out = zeros(Int, size(y_sat))
    if axis === :constant
        @inbounds for index in eachindex(out)
            value = y_sat[index]
            out[index] = (isnan(value) || value < threshold) ? 0 : 1
        end
        return out
    end
    isempty(sources) && throw(ArgumentError("$axis needs the source products to band on"))
    for source in sources
        size(source) == size(y_sat) ||
            throw(DimensionMismatch("every source product must have the anchor's shape"))
    end
    n_envelope = length(BLEND_ENVELOPE_EDGES)
    @inbounds for index in eachindex(out)
        agreement = 0
        envelope = -Inf
        for source in sources
            value = source[index]
            isnan(value) && continue
            envelope = max(envelope, value)
            value >= threshold && (agreement += 1)
        end
        agreement == 0 && continue
        agreement = min(agreement, 3)
        band = 1
        for (position, edge) in enumerate(BLEND_ENVELOPE_EDGES)
            envelope >= edge && (band = position)
        end
        out[index] = (agreement - 1) * n_envelope + band
    end
    return out
end

"""
What each cell of a banded blend leans toward. Every axis leans toward `fallback`
(`BLEND_FALLBACK_METHOD`), except that `:agreement_envelope_dry_raw` leans its dry band (the last
ID, no source product wet) toward the target's own anchor `y_sat`. A NaN anchor there leaves the
cell unblended, by `banded_blend_prediction`'s rule for a NaN fallback.

Returns `fallback` itself for every other axis, so those blends are unchanged.
"""
function blend_fallback(
    axis::Symbol, band::AbstractMatrix{Int}, fallback::AbstractMatrix{Float64},
    y_sat::AbstractMatrix{Float64},
)
    axis === :agreement_envelope_dry_raw || return fallback
    size(band) == size(fallback) == size(y_sat) || throw(DimensionMismatch(
        "band, fallback and y_sat must have the same shape",
    ))
    return ifelse.(band .== blend_band_count(axis), y_sat, fallback)
end

"""
Choose one blending weight per band on the inner selection split.

Scores the cells where the source and the fallback are both finite, as the scalar version does, so
the choice is not partly a question of which one dropped the harder cells. Each band's weight is
then minimised over `lambdas` independently, which is exact rather than a shortcut: see
`banded_blend_prediction` for why the bands separate.

Ties go to the smaller weight, matching `select_blend_lambda`: at equal error the answer that
departs less from the fitted model is the one to prefer, and it keeps the choice off the order of
the grid. A band with no scorable cell keeps weight 0, i.e. is left unblended - the same answer a
fold that cannot be helped gets.

Returns `nothing` when no cell is scorable at all, which the caller reports as a skipped fold.
"""
function select_blend_lambdas(
    y_obs::Matrix{Float64}, y_sat::Matrix{Float64}, band::Matrix{Int},
    prediction::Matrix{Float64}, fallback::Matrix{Float64}, lambdas::Vector{Float64},
    n_bands::Int; time_weights::Union{Nothing,Vector{Float64}}=nothing,
)
    shared = .!isnan.(y_obs) .& .!isnan.(y_sat) .&
        .!isnan.(prediction) .& .!isnan.(fallback)
    any(shared) || return nothing
    scored_band = ifelse.(shared, band, 0)
    chosen = zeros(Float64, n_bands)
    for target in 1:n_bands
        best_sse = Inf
        best_lambda = 0.0
        for lambda in lambdas
            sse = 0.0
            count = 0
            @inbounds for index in eachindex(scored_band)
                scored_band[index] == target || continue
                difference = (1 - lambda) * prediction[index] + lambda * fallback[index] -
                    y_obs[index]
                sse += difference * difference
                count += 1
            end
            count == 0 && continue
            if sse < best_sse - 1e-12
                best_sse = sse
                best_lambda = lambda
            end
        end
        chosen[target] = best_lambda
    end
    restricted_prediction = ifelse.(shared, prediction, NaN)
    restricted_fallback = ifelse.(shared, fallback, NaN)
    blended = banded_blend_prediction(
        scored_band, restricted_prediction, restricted_fallback, chosen,
    )
    best = _candidate_metrics(y_obs, y_sat, blended; time_weights)
    unblended = _candidate_metrics(y_obs, y_sat, restricted_prediction; time_weights)
    return (;
        lambda=chosen[1], lambdas=chosen, inner_RMSE=best.RMSE, inner_MAE=best.MAE,
        unblended_RMSE=unblended.RMSE, n=best.n, coverage=best.coverage,
        wet_cells=count(index -> scored_band[index] > 0, eachindex(scored_band)),
    )
end

"""
Choose the blending weight on the inner selection split.

Scores every candidate in `lambdas` on the cells where the source and the fallback are both finite,
so the choice is not partly a question of which one dropped the harder cells - the same discipline
`select_auto_method` applies to its contenders.

Ties go to the smaller `lambda`: at equal inner RMSE the answer that departs less from the fitted
model is the one to prefer, and it keeps the choice from depending on the order of the grid.

Returns `nothing` when no cell is scorable, which the caller reports as a skipped fold rather than
defaulting to a weight.
"""
function select_blend_lambda(
    y_obs::Matrix{Float64}, y_sat::Matrix{Float64}, prediction::Matrix{Float64},
    fallback::Matrix{Float64}, lambdas::Vector{Float64}, threshold::Float64;
    time_weights::Union{Nothing,Vector{Float64}}=nothing,
)
    shared = .!isnan.(y_obs) .& .!isnan.(y_sat) .&
        .!isnan.(prediction) .& .!isnan.(fallback)
    any(shared) || return nothing
    restricted_prediction = ifelse.(shared, prediction, NaN)
    restricted_fallback = ifelse.(shared, fallback, NaN)
    scored = NamedTuple[]
    for lambda in lambdas
        blended = satellite_wet_blend_prediction(
            y_sat, restricted_prediction, restricted_fallback, lambda, threshold,
        )
        push!(scored, merge((; lambda),
            _candidate_metrics(y_obs, y_sat, blended; time_weights)))
    end
    order = sortperm(scored; by=row -> (row.RMSE, row.MAE, row.lambda))
    best = scored[order[1]]
    unblended = scored[findfirst(row -> row.lambda == 0.0, scored)]
    return (;
        lambda=best.lambda, inner_RMSE=best.RMSE, inner_MAE=best.MAE,
        unblended_RMSE=unblended.RMSE, n=best.n, coverage=best.coverage,
        wet_cells=count(index -> shared[index] && y_sat[index] >= threshold,
            eachindex(shared)),
    )
end

"""
`matrix` moved one hour along `times`: `shift = -1` gives each cell its value at t-1, `+1` at t+1.

A cell whose neighbouring column is not exactly `shift` hours away (a gap in the common grid), or is
NaN there, keeps its own value, so a lagged stack feature degrades to the member it lags.
"""
function stack_hour_shift(matrix::AbstractMatrix{Float64}, times::AbstractVector, shift::Int)
    length(times) == size(matrix, 2) ||
        throw(DimensionMismatch("times must match the matrix's columns"))
    out = Matrix{Float64}(matrix)
    @inbounds for time in axes(matrix, 2)
        other = time + shift
        (1 <= other <= size(matrix, 2) && times[other] - times[time] == Hour(shift)) || continue
        for station in axes(matrix, 1)
            value = matrix[station, other]
            isnan(value) || (out[station, time] = value)
        end
    end
    return out
end

"""The stack's design, one matrix per entry of `STACK_FEATURES` and in that order."""
function stack_feature_matrices(
    mgwr::Matrix{Float64}, adw::Matrix{Float64}, raw::Matrix{Float64}, tps::Matrix{Float64},
    times::AbstractVector,
)
    size(mgwr) == size(adw) == size(raw) == size(tps) ||
        throw(DimensionMismatch("every stack member must have the same shape"))
    return [mgwr, adw, raw, tps,
        stack_hour_shift(mgwr, times, -1), stack_hour_shift(mgwr, times, 1),
        stack_hour_shift(adw, times, -1), stack_hour_shift(adw, times, 1)]
end

"""
Stack band per cell: `agreement + 1 + n_agreement * (distance_band - 1)`, in `1:stack_band_count()`.

`agreement` is `blend_band_matrix(:agreement_envelope, ...)` for these stations, whose band 0 (no
product wet) the stack fits like any other. `distance_km` is each station's distance to the nearest
station the predicting model was trained on, which is knowable at prediction time.
"""
function stack_band_matrix(agreement::AbstractMatrix{Int}, distance_km::AbstractVector{<:Real})
    length(distance_km) == size(agreement, 1) ||
        throw(DimensionMismatch("need one distance per station row"))
    n_agreement = blend_band_count(:agreement_envelope) + 1
    distance_band = [1 + count(edge -> distance >= edge, STACK_DISTANCE_EDGES)
                     for distance in distance_km]
    return agreement .+ 1 .+ n_agreement .* (distance_band .- 1)
end

"""
Non-negative least squares from normal equations, by enumerating the support.

Exact for the eight or so features the stack has (255 subsets of 8×8 solves). A singular subset is
skipped; if every subset is singular the zero vector comes back, which the caller never uses as
the fitted answer because the pooled system over every band is not singular in practice.
"""
function _nnls_normal(gram::Matrix{Float64}, rhs::Vector{Float64}, yy::Float64)
    p = length(rhs)
    best_sse, best = Inf, zeros(Float64, p)
    for support in 1:(2^p - 1)
        active = [j for j in 1:p if (support >> (j - 1)) & 1 == 1]
        solution = try
            cholesky(Symmetric(gram[active, active])) \ rhs[active]
        catch
            continue
        end
        all(>=(0), solution) || continue
        weights = zeros(Float64, p)
        weights[active] = solution
        sse = yy - 2 * dot(weights, rhs) + dot(weights, gram * weights)
        sse < best_sse && ((best_sse, best) = (sse, weights))
    end
    return best
end

"""
Fit one non-negative weight vector per stack band on the inner selection split.

Scores the cells where the gauge and every feature are finite, so no member's gaps decide which
cells a weight is fitted on. A band with fewer than `min_cells` such cells takes the pooled fit
over every band and is flagged in `fell_back`.

Returns `nothing` when no cell is usable, which the caller reports as a skipped fold. Otherwise
`weights` is `length(features) x n_bands`, beside the inner RMSE of the stack and of `features[1]`
alone on the same cells.
"""
function select_stack_weights(
    y_obs::Matrix{Float64}, features::Vector{Matrix{Float64}}, band::Matrix{Int}, n_bands::Int;
    min_cells::Int=STACK_MIN_BAND_CELLS,
)
    for feature in features
        size(feature) == size(y_obs) == size(band) || throw(DimensionMismatch(
            "every feature and the band matrix must have the observation's shape",
        ))
    end
    p = length(features)
    gram = zeros(Float64, p, p, n_bands)
    rhs = zeros(Float64, p, n_bands)
    yy = zeros(Float64, n_bands)
    used = zeros(Int, n_bands)
    row = zeros(Float64, p)
    @inbounds for index in eachindex(y_obs)
        observed = y_obs[index]
        isnan(observed) && continue
        usable = true
        for (j, feature) in enumerate(features)
            value = feature[index]
            isnan(value) && (usable = false; break)
            row[j] = value
        end
        usable || continue
        target = band[index]
        1 <= target <= n_bands || throw(ArgumentError("band $target is outside 1:$n_bands"))
        used[target] += 1
        yy[target] += observed * observed
        for i in 1:p
            rhs[i, target] += row[i] * observed
            for j in i:p
                gram[i, j, target] += row[i] * row[j]
            end
        end
    end
    sum(used) == 0 && return nothing
    for target in 1:n_bands, i in 1:p, j in 1:(i - 1)
        gram[i, j, target] = gram[j, i, target]
    end
    pooled = _nnls_normal(dropdims(sum(gram; dims=3); dims=3), vec(sum(rhs; dims=2)), sum(yy))
    weights = zeros(Float64, p, n_bands)
    fell_back = falses(n_bands)
    sse_stack, sse_source = 0.0, 0.0
    for target in 1:n_bands
        fell_back[target] = used[target] < min_cells
        weights[:, target] = fell_back[target] ? pooled :
            _nnls_normal(gram[:, :, target], rhs[:, target], yy[target])
        used[target] == 0 && continue
        w = weights[:, target]
        sse_stack += yy[target] - 2 * dot(w, rhs[:, target]) + dot(w, gram[:, :, target] * w)
        sse_source += yy[target] - 2 * rhs[1, target] + gram[1, 1, target]
    end
    n = sum(used)
    return (; weights, used, fell_back, n,
        inner_RMSE=sqrt(max(sse_stack, 0.0) / n), source_RMSE=sqrt(max(sse_source, 0.0) / n))
end

"""
Apply per-band stack weights. A cell with any feature missing keeps `features[1]`'s value - the
source prediction, NaN if the source itself is - so a gap in a secondary member can neither shrink
the stack's coverage below its source's nor poison a cell, the rule the blends follow too.
"""
function stacked_prediction(
    features::Vector{Matrix{Float64}}, band::Matrix{Int}, weights::Matrix{Float64},
)
    size(weights, 1) == length(features) ||
        throw(DimensionMismatch("need one weight per feature"))
    source = first(features)
    out = similar(source)
    @inbounds for index in eachindex(source)
        total = 0.0
        complete = true
        for (j, feature) in enumerate(features)
            value = feature[index]
            isnan(value) && (complete = false; break)
            total += weights[j, band[index]] * value
        end
        out[index] = complete ? total : source[index]
    end
    return out
end
