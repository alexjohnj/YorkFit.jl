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
function fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σXs::Vector{T}, σYs::Vector{T}, rs::Vector{T}; iter=20, tol=1E-6)
    println(iter, tol)
end

end # module
