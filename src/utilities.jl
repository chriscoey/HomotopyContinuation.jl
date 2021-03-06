filter_kwargs(predicate, kwargs) = filter( x -> predicate(first(x)), kwargs)

# we define our own widen
_widen(T) = widen(T)
# Wait until DoubleFloat64 is released
# _widen(::Type{Float64}) = FastDouble

affine(xs::AbstractVector) = xs[2:end] ./ x[1]

"""
    embed_projective_if_necessary(x, H)

Embeds a vector into the projective space if necessary, i.e. if it's length is one less
than the number of variables of `H`. `H` is assumed to be homogenized. After the (eventual)
embedding the value is normalized.
"""
function embed_projective_if_necessary!(x::AbstractVector{T}, H::AbstractHomotopy{T}) where T
    N = Homotopies.nvariables(H)
    n = length(x)
    if N - 1 == n
        unshift!(x, one(T))
    elseif N != n
        error("A start value has length $n. Excepted length $N or $(N-1).")
    end
    x
end

"""
    projectivenorm2(a, b)

Calculate the squared infinity norm of |a-b|, but first brings `a` and `b` on the same patch.
Brings fist `a` and `b` on the same patch by finding diving `a` through it's maximum value
(w.r.t. to the absolute value) with index `i` and diving `b` through `b[i]`.
Then computes the norm of the differences.
"""
function projectivenorm2(a::AbstractVector{T}, b::AbstractVector{T}) where T
    maxind = 1
    maxval = abs2(first(a))
    for i = 2:length(a)
        val = abs2(a[i])
        if val > maxval
            maxind = i
            maxval = val
        end
    end
    out = real(zero(T))
    adiv = a[maxind]
    bdiv = b[maxind]

    for i=1:length(a)
        out = max(out, abs2(a[i] / adiv - b[i] / bdiv))
    end
    out
end


"""
    UnitRootsIterator(r, n)

Construct an infinite iterator which returns the `n`-th scaled unit roots, i.e.
the values ``r⋅exp(i2πk/n)`` for ``k=0,1,...``.
"""
struct UnitRootsIterator
    radius::Float64
    order::Float64
end
UnitRootsIterator(r::Real, order::Real) = UnitRootsIterator(float(r), float(order))

Base.start(::UnitRootsIterator) = 0
Base.next(loop::UnitRootsIterator, k::Int) = (loop.radius *  exp(im * 2π * k / loop.order), k + 1)
Base.done(::UnitRootsIterator, ::Int) = false
Base.iteratorsize(::UnitRootsIterator) = Base.IsInfinite()
Base.iteratoreltype(::UnitRootsIterator) = Base.HasEltype()
Base.eltype(::UnitRootsIterator) = Complex128
