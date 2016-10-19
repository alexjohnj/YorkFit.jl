module YorkFit

# package code goes here
"""Calculate `a` and `b` using the method of least squares. Returns the tuple
`(a, b)` where `a` is the intercept and `b` is the gradient."""
function lsq{T<:Real}(Xs::Vector{T}, Ys::Vector{T})::Tuple{T,T}
    A = hcat(ones(length(Xs)), Xs)
    m = (A' * A) \ A' * Ys

    (m[1], m[2])
end

end # module
