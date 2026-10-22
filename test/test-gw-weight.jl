using Test, MixedGWR
using MixedGWR: gw_weight_bisquare

@testset "gw_weight: one adaptive and fixed rule" begin
    bisquare(d, b) = d < b ? (1 - (d / b)^2)^2 : 0.0
    distances = [5.0, 1.0, 9.0, 3.0, Inf, 7.0]
    # Adaptive: the bandwidth is the k-th smallest finite distance (here 5.0).
    weights = gw_weight(distances, 3.0; kernel=BISQUARE)
    @test weights == [isfinite(d) ? bisquare(d, 5.0) : 0.0 for d in distances]
    @test count(>(0), weights) == 2
    # A kernel given as a function weighs exactly like its index.
    @test gw_weight(distances, 3.0; kernel=gw_weight_bisquare) == weights
    # More neighbours than finite distances: GWmodel's k/m * maximum over the finite ones.
    @test gw_weight(distances, 10.0; kernel=BISQUARE) ==
        [isfinite(d) ? bisquare(d, 10 / 5 * 9.0) : 0.0 for d in distances]
    # A non-finite distance never gets weight, even at an infinite fixed bandwidth.
    @test gw_weight(distances, Inf; kernel=BISQUARE, adaptive=false) == Float64.(isfinite.(distances))
    @test all(iszero, gw_weight(fill(Inf, 4), 2.0; kernel=BISQUARE))
    # A zero bandwidth (k coincident points) degenerates to exact coincidence, for any kernel.
    @test gw_weight([0.0, 0.0, Inf], 2.0; kernel=BISQUARE) == [1.0, 1.0, 0.0]
    @test gw_weight([0.0, 0.0, 3.0], 2.0; kernel=GAUSSIAN) == [1.0, 1.0, 0.0]
    @test gw_weight(zeros(4), 2.0; kernel=BISQUARE) == ones(4)
    # The matrix form is the vector form per column.
    D = [distances reverse(distances)]
    @test gw_weight(D, 3.0; kernel=BISQUARE) ==
        hcat(gw_weight(D[:, 1], 3.0; kernel=BISQUARE), gw_weight(D[:, 2], 3.0; kernel=BISQUARE))
end

@testset "haversine_distance_matrix" begin
    lonlat = [110.0 32.0; 110.5 32.0; 111.0 33.0]
    D = haversine_distance_matrix(lonlat, lonlat)
    @test all(iszero, D[i, i] for i in 1:3)
    @test D ≈ D'
    # One degree of latitude on the mean-radius sphere.
    @test haversine_distance_matrix([110.0 32.0], [110.0 33.0])[1, 1] ≈ MixedGWR.EARTH_RADIUS_KM * pi / 180
    L = loocv_distances(D)
    @test all(isinf(L[i, i]) for i in 1:3)
    @test all(L[i, j] == D[i, j] for i in 1:3, j in 1:3 if i != j)
end
