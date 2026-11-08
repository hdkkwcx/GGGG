# ------------------------------------------- D7: one sample for every product (matched cells)

"""
Methods whose out-of-fold value must be finite in a matched cell: the mask-defining methods, plus
`auto` and the blends that a paper table reports beside them. A method absent from a product's
directory is skipped, so a run without blends needs no shorter list.
"""
const MATCHED_METHODS = vcat(MASK_METHODS, [
    "auto", "blend_residual_gwr", "blend_mixed_gwr", "blend_mgwr",
    "blend_agrenv_residual_gwr", "blend_agrenv_mixed_gwr", "blend_agrenv_mgwr",
])

"""
The cells every product of a run is scored on, in one scheme: `(; matched, masks)`.

`masks[p]` is product `p`'s own `common_evaluation_mask.csv`, which already needs the gauge and
every `MASK_METHODS` prediction finite (`raw` among them, so the product's own value). `matched` is
their intersection over `products`, also requiring every other method in `methods` to be finite for
every product. Each product is otherwise scored on its own mask, so "which product is better" would
compare different station-hours; on `matched` every product and method sees the same ones.

Files are read one at a time, each checked hour by hour against `times`
(`DayBlockBootstrap.read_wide_mask` / `read_wide_values`).
"""
function matched_evaluation_mask(scheme_dir::AbstractString, products::Vector{String},
        ids::Vector{String}, times::Vector{DateTime}; methods::Vector{String}=MATCHED_METHODS)
    isempty(products) && throw(ArgumentError("no products to match"))
    matched = trues(length(ids), length(times))
    masks = Dict{String,BitMatrix}()
    for product in products
        product_dir = joinpath(scheme_dir, lowercase(product))
        masks[product] = read_wide_mask(
            joinpath(product_dir, "common_evaluation_mask.csv"), ids, times)
        matched .&= masks[product]
        for method in setdiff(methods, MASK_METHODS)
            path = joinpath(product_dir, "oof_$(method).csv")
            isfile(path) || continue
            matched .&= .!isnan.(read_wide_values(path, ids, times))
        end
    end
    return (; matched=BitMatrix(matched), masks)
end

"""
The `run_status.csv` rows that would take a whole fold out of the matched sample: a method in
`methods` that `failed` (or was `skipped` because the method it builds on failed) for one of
`products` in `scheme`. In the matched sample such a fold is lost for every product, so callers stop
on a non-empty result instead of shrinking the sample silently. `partial` folds are allowed: their
missing cells drop out cell by cell, and `matched_sample_table` reports what they cost.
"""
matched_fold_failures(status::DataFrame; scheme::AbstractString, products::Vector{String},
                      methods::Vector{String}=MATCHED_METHODS) =
    filter(row -> row.scheme == scheme && row.product in products && row.method in methods &&
                  row.status in ("failed", "skipped"), status)

"""
Per product and fold (and `"all"`), the cells in the product's own mask, in the matched sample, and
the share kept, with the hours holding at least one matched cell.
"""
function matched_sample_table(masks::AbstractDict{String,BitMatrix}, matched::BitMatrix,
                              fold_of::Vector{Int}; scheme::AbstractString)
    rows = NamedTuple[]
    folds = vcat([("all", collect(eachindex(fold_of)))],
        [(string(fold), findall(==(fold), fold_of)) for fold in sort(unique(fold_of))])
    for product in sort(collect(keys(masks))), (fold, stations) in folds
        n_own = count(view(masks[product], stations, :))
        n_matched = count(view(matched, stations, :))
        push!(rows, (; scheme=String(scheme), product, fold, n_own, n_matched,
            matched_share_of_own=n_own == 0 ? NaN : n_matched / n_own,
            hours_matched=count(vec(any(view(matched, stations, :); dims=1)))))
    end
    return DataFrame(rows)
end
