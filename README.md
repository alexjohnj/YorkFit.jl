# YorkFit.jl

[![Build Status](https://travis-ci.org/alexjohnj/YorkFit.jl.svg?branch=master)](https://travis-ci.org/alexjohnj/YorkFit.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/wbd43g85rvkk56sc?svg=true)](https://ci.appveyor.com/project/alexjohnj/yorkfit-jl)


This package implements the line fitting algorithm described by York et
al. (2004). It provides a method for calculating the line that fits data with
errors in the x and y variables. YorkFit exposes a single function

``` julia
fit{T<:Real}(Xs::Vector{T}, Ys::Vector{T}, σXs::Vector{T}, σYs::Vector{T}, rs::Vector{T}; kwargs...)::Tuple{T,T,T,T}
```

This fits a line to the data in the vectors `Xs` and `Ys` with the errors `𝞂Xs`
and `𝞂Ys` that have the correlation coefficients in `rs`. The function returns a
4-tuple containing the y-intercept, gradient, error in the intercept and error
in the gradient. Several methods for `fit` are available that, for example, only
require a single `𝞂X` value. The function is documented so entering `?fit` in
the REPL will give more information on the available methods.

## Installation

This isn't available in Julia's package repository yet. For now run

``` julia
Pkg.clone("https://travis-ci.org/alexjohnj/YorkFit.jl.git")
```

to get YorkFit.

## References

York, Derek, et al. "Unified equations for the slope, intercept, and standard
errors of the best straight line." _American Journal of Physics_ 72.3 (2004):
367-375.
