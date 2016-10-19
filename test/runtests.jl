using YorkFit
using Base.Test

# Test least-squares estimate
let
    Xs = collect(0:10.0)
    Ys = 2.0 * Xs + 3.0
    (a, b) = YorkFit.lsq(Xs, Ys)
    @test isapprox(a, 3.0)
    @test isapprox(b, 2.0)
end

