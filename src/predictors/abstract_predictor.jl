"""
    AbstractPredictor

The path tracking problem can be considered a initial value problem.
A predictor make use of this fact and predicts for a given pair ``(x₀,t)``
with ``H(x₀,t)=0`` a new ``x`` such that ``H(x, t + Δt)=0``.

The differential equation is
```math
x′(t) = -Hₓ(x(t), t)⁻¹Hₜ(x(t), t)
```
which follows from the fact ``∂tH(x(t),t) ≡ 0 ∀ t∈[0,1]`` and the total derivative
of ``H`` w.r.t. ``t``.
"""
abstract type AbstractPredictor end
abstract type AbstractPredictorCache end

"""
    cache(::AbstractPredictor, ::HomotopyWithCache, x, t)::AbstractPredictorCache

Construct a cache to avoid allocations.
"""
function cache end


"""
    predict!(xnext, ::AbstractPredictor, ::AbstractPredictorCache, H::HomotopyWithCache, x, t, Δt)

Perform a prediction step for the value of `x` with step size `Δt`.
"""
function predict! end

struct NullPredictor <: AbstractPredictor end
struct NullPredictorCache <: AbstractPredictorCache end

cache(::NullPredictor, H, x, t) = NullPredictorCache()

function predict!(xnext, ::NullPredictor, ::NullPredictorCache, H, x, t, dt)
    xnext .= x
    nothing
end
