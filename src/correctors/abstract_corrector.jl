# Interface
abstract type AbstractCorrector end
abstract type AbstractCorrectorCache end

"""
    Result{T}

Structure holding information about a `correct!` step. The fields are
* `converged::Bool` Indicating whether the correction was successfull
* `res::T` The residual of the final result.
* `iters::Int` The number of iterations used.
"""
struct Result{T}
    converged::Bool
    res::T
    iters::Int
end

"""
    cache(::AbstractCorrector, ::HomotopyWithCache{M, N}, x, t)::AbstractCorrectorCache

Construct a cache to avoid allocations.
"""
function cache end


"""
    correct!(xnext, ::AbstractCorrector, ::AbstractCorrectorCache, H::HomotopyWithCache, x, t, tol)::Result

Perform a correction step such that in the end `H(xnext,t) < tol`.
Returns a [`Result`](@ref).
"""
function correct! end
