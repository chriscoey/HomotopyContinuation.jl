struct Euler <: AbstractPredictor end
struct EulerCache{T} <: AbstractPredictorCache
    A::Matrix{T}
    b::Vector{T}
end

cache(::Euler, H, x, t) = EulerCache(jacobian(H, x, t), dt(H, x, t))

"""
    minus_x_prime!(out, H, x, t, A)

Evaluate `Hₓ(x(t), t)⁻¹Hₜ(x(t), t)` and store the result in `out`. `A` needs
to be able to store the Jacobian of `H`.
"""
function minus_x_prime!(out, H, x, t, A)
    jacobian_and_dt!(A, out, H, x, t)
    solve_with_lu_inplace!(A, out)
    out
end


function predict!(xnext, ::Euler, cache::EulerCache, H::HomotopyWithCache, x, t, Δt)
    minus_x_prime!(cache.b, H, x, t, cache.A)
    @. xnext = x - Δt * cache.b
    nothing
end
