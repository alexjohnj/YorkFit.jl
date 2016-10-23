module YorkFit

"""Calculate the weight of the error `σX`."""
ω(σX::Real)::Real = σX^-2
α{T<:Real}(ωX::T, ωY::T)::T = sqrt(ωX * ωY)
W{T<:Real}(ωX::T, ωY::T, r::T, b::T, α::T)::T =
    (ωX * ωY) / (ωX + b^2 * ωY - 2 * b * r * α)

"""Compute Σ(Ws * Xs) / Σ(Xs)"""
barvar{T<:Real}(Ws::Vector{T}, Xs::Vector{T})::T = sum(Ws .* Xs) / sum(Ws)

function β{T<:Real}(W::T, U::T, V::T, ωX::T, ωY::T, b::T, r::T, α)::T
    W * ((U / ωY) + (b * V / ωX) - (r / α) * (b * U + V))
end

"""Calculate `a` and `b` using the method of least squares. Returns the tuple
`(a, b)` where `a` is the intercept and `b` is the gradient."""
function lsq{T<:Real}(Xs::Vector{T}, Ys::Vector{T})::Tuple{T,T}
    A = hcat(ones(length(Xs)), Xs)
    m = (A' * A) \ A' * Ys

    (m[1], m[2])
end

"""
    fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σX::T, σYs::Vector{T}, rs::Vector{T}; kwargs...)
    fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σX::T, σYs::Vector{T}, r::T=T(0); kwargs...)
    fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σXs::Vector{T}, σY::T, rs::Vector{T}; kwargs...)
    fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σXs::Vector{T}, σY::T, r::T=T(0); kwargs...)
    fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σX::T, σY::T, rs::Vector{T}; kwargs...)
    fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σX::T, σY::T, r::T=T(0); kwargs...)
    fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σXs::Vector{T}, σYs::Vector{T}, r::T=T(0); kwargs...)

    fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σXs::Vector{T}, σYs::Vector{T}, rs::Vector{T}; kwargs...)::Tuple{T,T,T,T}

Fits a slope and intercept to a set of observations using the method of York et
al. (2004).

Arguments
=========

- `Xs` -- A vector of observations.
- `Ys` -- A vector of observations related to `Xs` by `Ys = a + b * Xs`.
- `σXs` -- Errors in the observations in `Xs`. A vector of equal length to `Xs` or
  a scalar value if all the errors are equal.
- `σYs` -- Errors in the observations in `Ys`. A vector of equal length to `Ys` or
  a scalar value if all the errors are equal.
- `rs` -- Correlation coefficient between errors. A vector of equal length to
  `Xs` and `Ys` or a scalar value if the coefficient is constant. Defaults to 0.

Keyword Arguments
==================

- `niter=400` -- Number of iterations to calculate the slope `b` with.
- `tol=1E-15` -- Stopping criteria. When the difference between two successive
  iterations of `b` falls below this, the iterations stop.

Returns
=======

A 4-tuple containing `(a, b, σa, σb)` where `a` is the intercept, `b` is the
slope and `σa` and `σb` are the errors in the intercept and slope.

References
==========

York, D., Evensen, N.M., Martınez, M.L. and Delgado, J.D.B., 2004. Unified
equations for the slope, intercept, and standard errors of the best straight
line. American Journal of Physics, 72(3), pp.367-375.
"""
fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σX::T, σYs::Vector{T}, rs::Vector{T}; kwargs...) =
    fit(Xs, Ys, fill(σX, length(Xs)), σYs, rs; kwargs...)
fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σX::T, σYs::Vector{T}, r::T=T(0); kwargs...) =
    fit(Xs, Ys, fill(σX, length(Xs)), σYs, fill(r, length(Xs)); kwargs...)
fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σXs::Vector{T}, σY::T, rs::Vector{T}; kwargs...) =
    fit(Xs, Ys, σXs, fill(σY, length(Ys)), rs; kwargs...)
fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σXs::Vector{T}, σY::T, r::T=T(0); kwargs...) =
    fit(Xs, Ys, σXs, fill(σY, length(Ys)), fill(r, length(Xs)); kwargs...)
fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σX::T, σY::T, rs::Vector{T}; kwargs...) =
    fit(Xs, Ys, fill(σX, length(Xs)), fill(σY, length(Ys)), rs; kwargs...)
fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σX::T, σY::T, r::T=T(0); kwargs...) =
    fit(Xs, Ys, fill(σX, length(Xs)), fill(σY, length(Ys)), fill(r, length(Xs)); kwargs...)
fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σXs::Vector{T}, σYs::Vector{T}, r::T=T(0); kwargs...) =
    fit(Xs, Ys, σXs, σYs, fill(r, length(Xs)); kwargs...)
function fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σXs::Vector{T}, σYs::Vector{T}, rs::Vector{T}; niter=400, tol=1E-15)::Tuple{T,T,T,T}
    # TODO: Check for dimension mismatch of vectors
    (a, b) = lsq(Xs, Ys)

    ωXs = ω.(σXs)
    ωYs = ω.(σYs)
    αs = α.(ωXs, ωYs)

    Ws = zeros(length(Xs))
    βs = zeros(length(Xs))
    barX = 0
    barY = 0

    for i in 1:niter
        Ws = W.(ωXs, ωYs, rs, b, αs)
        barX = barvar(Ws, Xs)
        barY = barvar(Ws, Ys)
        Us = Xs - barX
        Vs = Ys - barY
        βs = β.(Ws, Us, Vs, ωXs, ωYs, b, rs, αs)

        b_old = b
        b = sum(Ws .* βs .* Vs) / sum(Ws .* βs .* Us)

        if abs((b - b_old) / b) < tol
            break
        end
    end

    a = barY - b * barX
    # Adjusted x values
    xs = barX + βs
    barx = sum(Ws .* xs) / sum(Ws)
    us = xs - barx

    σb = sqrt(sum(Ws .* us.^2)^-1)
    σa = sqrt(sum(Ws)^-1 + barx^2 * σb^2)

    return (a, b, σa, σb)
end

end # module
