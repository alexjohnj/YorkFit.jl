using YorkFit
using Base.Test

# Test least-squares estimate
let
    Xs = collect(1:10.0)
    Ys = 2.0 * Xs + 3.0
    (a, b) = YorkFit.lsq(Xs, Ys)
    @test isapprox(a, 3.0)
    @test isapprox(b, 2.0)
end

# Test components of York method with simple cases
let
    Xs = collect(1:10.0)
    Ys = 2.0*Xs + 3.0
    σXs = fill(0.05, 10)
    σYs = fill(0.02, 10)
    r = 0.0

    ωXs = YorkFit.ω.(σXs)
    ωYs = YorkFit.ω.(σYs)

    @test isapprox(ωXs, fill(400.0, 10))
    @test isapprox(ωYs, fill(2500.0, 10))

    αs = YorkFit.α.(ωXs, ωYs)
    @test isapprox(αs, fill(1000.0, 10))

    Ws = YorkFit.W.(ωXs, ωYs, fill(r, 10), 2.0, αs)
    @test isapprox(Ws, fill(1250/13, 10))

    barX = YorkFit.barvar(Ws, Xs)
    barY = YorkFit.barvar(Ws, Ys)

    @test isapprox(barX, 5.5)
    @test isapprox(barY, 14.0)

    Us = Xs - barX
    Vs = Ys - barY
    βs = YorkFit.β.(Ws, Us, Vs, ωXs, ωYs, 2.0, r, αs)
    @test isapprox(βs, collect(-4.5:1:4.5))
end
