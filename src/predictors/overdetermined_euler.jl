
struct OverdeterminedEuler <: AbstractPredictor end

struct OverdeterminedEulerCache{ aType<:AbstractMatrix, asub, bType<:AbstractVector, bsub, cType<:AbstractVector} <: AbstractPredictorCache
    A::aType
    A_sub::asub # this is a MxN matrix
    b::bType
    b_sub::bsub
    c::cType# this is a N-Vector
    m::Int
end

function cache(::OverdeterminedEuler, H, x, t)
    J = jacobian(H, x, t)
    n, m = size(J)
    A = zeros(eltype(J), n+1, m)
    A_sub = @view A[1:n, :]
    b = zeros(eltype(J), n+1)
    b_sub = @view b[1:n]
    c = zeros(eltype(J), m)
    OverdeterminedEulerCache{typeof(A), typeof(A_sub), typeof(b), typeof(b_sub), typeof(c)}(A, A_sub, b, b_sub, c, m)
end

function predict!(xnext, ::OverdeterminedEuler, cache::OverdeterminedEulerCache, H::HomotopyWithCache, x, t, Δt)

    jacobian_and_dt!(cache.A_sub, cache.b_sub, H, x, t)
    for j=1:cache.m
        cache.A[end, j] = conj(x[j])
    end
    cache.b[end] = 0

    A_ldiv_B!(qrfact!(cache.A), cache.b)
    for i=1:cache.m
        cache.c[i] = cache.b[i]
    end
    @. xnext = x - Δt * cache.c
    nothing
end
