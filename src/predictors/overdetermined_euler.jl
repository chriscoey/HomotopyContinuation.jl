
struct OverdeterminedEuler <: AbstractPredictor end

struct OverdeterminedEulerCache{ aType<:AbstractMatrix, asub, bType<:AbstractVector, bsub} <: AbstractPredictorCache
    A::aType
    A_sub::asub # this is a MxN matrix
    b::bType
    b_sub::bsub # this is a N-Vector
    m::Int
end

function cache(::OverdeterminedEuler, H, x, t)
    n = length(H)
    m = nvariables(H)
    A = zeros(eltype(A), n+1, m)
    A_sub = @view A[1:n, :]
    b = zeros(eltype(x), n+1)
    b_sub = @view b[1:n]
    OverdeterminedEulerCache{typeof(A), typeof(A_sub), typeof(b), typeof(b_sub)}(A, A_sub, b, b_sub, m)
end

function predict!(xnext, ::OverdeterminedEuler, cache::OverdeterminedEulerCache, H::HomotopyWithCache, x, t, Δt)

    jacobian_and_dt!(cache.A_sub, cache.b_sub, H, x, t)
    for j=1:cache.m
        cache.A[end, j] = conj(x[j])
    end
    cache.b[end] = zero(T)

    A_ldiv_B!(qrfact!(cache.A), cache.b)
    @. xnext = x - Δt * cache.b
    nothing
end
