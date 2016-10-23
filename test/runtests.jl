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

# Test fit methods
let
    Xs = collect(1.0:10.0)
    Ys = 2.0*Xs + 3.0
    σX = 0.25
    σY = 0.5
    r = 0.25

    # Expected results
    Ea = 3.0
    Eb = 2.0
    Eσb = 0.06741998624
    Eσa = 0.41833001326

    function testvals(a, b, σa, σb)
        @test isapprox(a, Ea)
        @test isapprox(b, Eb)
        @test isapprox(σa, Eσa)
        @test isapprox(σb, Eσb)
    end

    # Vector r, constant σX, Vector σY
    (a, b, σa, σb) = YorkFit.fit(Xs, Ys, σX, fill(σY, 10), fill(r, 10), niter=2)
    testvals(a, b, σa, σb)

    # Constant r, constant σX, Vector σY
    (a, b, σa, σb) = YorkFit.fit(Xs, Ys, σX, fill(σY, 10), r, niter=2)
    testvals(a, b, σa, σb)

    # Vector r, vector σX, constant σY
    (a, b, σa, σb) = YorkFit.fit(Xs, Ys, fill(σX, 10), σY, fill(r, 10), niter=2)
    testvals(a, b, σa, σb)

    # Constant r, vector σX, constant σY
    (a, b, σa, σb) = YorkFit.fit(Xs, Ys, fill(σX, 10), σY, r, niter=2)
    testvals(a, b, σa, σb)

    # Vector r, constant σX, σY
    (a, b, σa, σb) = YorkFit.fit(Xs, Ys, σX, σY, fill(r, 10), niter=2)
    testvals(a, b, σa, σb)

    # Constant σX, σY, r
    (a, b, σa, σb) = YorkFit.fit(Xs, Ys, σX, σY, r, niter=2)
    testvals(a, b, σa, σb)

    # Constant r, vector σX, Vector σY
    (a, b, σa, σb) = YorkFit.fit(Xs, Ys, fill(σX, 10), fill(σY, 10), r, niter=2)
    testvals(a, b, σa, σb)

    # Vector r, vector σX, Vector σY
    (a, b, σa, σb) = YorkFit.fit(Xs, Ys, fill(σX, 10), fill(σY, 10), fill(r, 10), niter=2)
    testvals(a, b, σa, σb)
end
