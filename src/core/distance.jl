export haversine_distance_matrix, loocv_distances

# Every station distance in the project is on a sphere of this radius: the IUGG mean Earth
# radius, in km.
const EARTH_RADIUS_KM = 6371.0088

"""
Pairwise great-circle distances in km between `(lon, lat)` rows in degrees, with rows=train and
columns=target.
"""
function haversine_distance_matrix(
    train_lonlat::AbstractMatrix{<:Real}, target_lonlat::AbstractMatrix{<:Real},
)
    ntrain, ntarget = size(train_lonlat, 1), size(target_lonlat, 1)
    result = Matrix{Float64}(undef, ntrain, ntarget)
    # Each row's radians and cos(lat) depend on that row alone, so they are computed once.
    train_lon = [deg2rad(Float64(train_lonlat[i, 1])) for i in 1:ntrain]
    train_lat = [deg2rad(Float64(train_lonlat[i, 2])) for i in 1:ntrain]
    train_coslat = cos.(train_lat)
    @inbounds for j in 1:ntarget
        lon2, lat2 = deg2rad(Float64(target_lonlat[j, 1])), deg2rad(Float64(target_lonlat[j, 2]))
        coslat2 = cos(lat2)
        for i in 1:ntrain
            dlon, dlat = lon2 - train_lon[i], lat2 - train_lat[i]
            a = sin(dlat / 2)^2 + train_coslat[i] * coslat2 * sin(dlon / 2)^2
            result[i, j] = 2 * EARTH_RADIUS_KM * asin(sqrt(clamp(a, 0.0, 1.0)))
        end
    end
    return result
end

"""A copy of `distances` with `Inf` on the diagonal, which leaves each target out of its own fit."""
function loocv_distances(distances::AbstractMatrix{<:Real})
    d = Matrix{Float64}(distances)
    for i in 1:min(size(d)...)
        d[i, i] = Inf
    end
    return d
end
