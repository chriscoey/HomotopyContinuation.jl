
# Newton
"""
    OverdeterminedNewton()

A classical simple Newton operator for overdetermined systems using the LU factorization
to solve the linear systems.
"""
struct OverdeterminedNewton <: AbstractCorrector end

struct OverdeterminedNewtonCache{ aType<:AbstractMatrix, asub, bType<:AbstractVector, bsub} <: AbstractCorrectorCache
    A::aType
    A_sub::asub # this is a MxN matrix
    b::bType
    b_sub::bsub # this is a N-Vector
    m::Int
end

function cache(::OverdeterminedNewton, H::HomotopyWithCache, x, t)
    n = length(H)
    m = nvariables(H)
    A = zeros(eltype(A), n+1, m)
    A_sub = @view A[1:n, :]
    b = zeros(eltype(x), n+1)
    b_sub = @view b[1:n]
    OverdeterminedNewtonCache{typeof(A), typeof(A_sub), typeof(b), typeof(b_sub)}(A, A_sub, b, b_sub, m)
end

function correct!(xnext, ::OverdeterminedNewton,
    cache::OverdeterminedNewtonCache,
    H::HomotopyWithCache,
    x,
    t,
    tol,
    maxiters)

    A, A_sub, b, b_sub, m = cache.A, cache.A_sub, cache.b, cache.b_sub, cache.m

    k = 0
    # println("start")
    while true
        k += 1

        evaluate!(b_sub, H, xnext, t)
        b[end] = zero(T)

        if infinity_norm(b) < abstol
            return true
        elseif k > maxiters
            return false
        end

        # put jacobian in A
        jacobian!(A_sub, H, xnext, t)
        for j=1:m
            A[end, j] = conj(xnext[j])
        end

        A_ldiv_B!(qrfact!(A), b)

        @. xnext = xnext - b
        normalize!(xnext)
    end
end
