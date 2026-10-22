# GWR against R's GWmodel on the Shiyan station fixture. The only test file that needs R (with
# GWmodel installed); test-shiyan-fixture.jl includes it only when the fixture exists.
using RCall
R"""
library(GWmodel)
"""

@testset "GWR" begin
  for adaptive in [true, false]
    for kernel in 0:4
      r = R"GWmodel:::gwr_q($x1, $y, $dMat, 20.0, $kernel + 1, $adaptive)" |> rcopy
      jl = GWR(x1, y, dMat, 20.0; kernel, adaptive)
      @test r ≈ jl
    end
  end
end
