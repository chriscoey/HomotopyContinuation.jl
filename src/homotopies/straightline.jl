# StraightLineHomotopy implementation
"""
    StraightLineHomotopy(G, F; gamma=exp(i * 2π*rand()))

Construct the homotopy ``H(x, t) = γtG(x) + (1-t)F(x)``.
"""
struct StraightLineHomotopy{S<:AbstractSystem,T<:AbstractSystem} <: AbstractStartTargetHomotopy
    start::S
    target::T
    gamma::Complex{Float64}

end
function StraightLineHomotopy(start::AbstractSystem, target::AbstractSystem; gamma=randomish_gamma())
    StraightLineHomotopy(start, target, gamma)
end

function randomish_gamma()
    # Usually values near 1, i, -i, -1 are not good randomization
    # Therefore we artificially constrain the choices
    theta = rand() * 0.30 + 0.075 + (rand(Bool) ? 0.0 : 0.5)
    cis(2π * theta)
end

start(H::StraightLineHomotopy) = H.start
target(H::StraightLineHomotopy) = H.target

Base.size(H::StraightLineHomotopy) = size(H.start)

"""
    gamma(H::StraightLineHomotopy)

Obtain the gamma used in the StraightLineHomotopy.
"""
Base.Math.gamma(H::StraightLineHomotopy) = H.gamma

"""
    γ(H)

Obtain the gamma used in the StraightLineHomotopy.
"""
γ(H::StraightLineHomotopy) = gamma(H)

function evaluate!(u, H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    Systems.evaluate!(c.u, start(H), x, c.start)
    Systems.evaluate!(u, target(H), x, c.target)

    u .= (γ(H) * t) .* c.u .+ (1 - t) .* u

    u
end
function evaluate(H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    G = Systems.evaluate(start(H), x, c.start)
    F = Systems.evaluate(target(H), x, c.target)
    (γ(H) * t) * G + (1 - t) * F
end
(H::StraightLineHomotopy)(x, t, c=cache(H, x, t)) = evaluate(H, x, t, c)

function dt!(u, H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    Systems.evaluate!(c.u, start(H), x, c.start)
    Systems.evaluate!(u, target(H), x, c.target)

    u .= γ(H) .* c.u .- u

    u
end
function dt(H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    G = Systems.evaluate(start(H), x, c.start)
    F = Systems.evaluate(target(H), x, c.target)
    γ(H) .* G .- F
end

function jacobian!(U, H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    Systems.jacobian!(c.U, start(H), x, c.start)
    Systems.jacobian!(U, target(H), x, c.target)

    U .= (γ(H) * t) .* c.U .+ (1 - t) .* U

    U
end
function jacobian(H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    G = Systems.jacobian(start(H), x, c.start)
    F = Systems.jacobian(target(H), x, c.target)
    (γ(H) * t) .* G .+ (1 - t) .* F
end

function evaluate_and_jacobian!(u, U, H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    Systems.evaluate_and_jacobian!(c.u, c.U, start(H), x, c.start)
    Systems.evaluate_and_jacobian!(u, U, target(H), x, c.target)

    u .= (γ(H) * t) .* c.u .+ (1 - t) .* u
    U .= (γ(H) * t) .* c.U .+ (1 - t) .* U

    nothing
end

function jacobian_and_dt!(U, u, H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    Systems.evaluate_and_jacobian!(c.u, c.U, start(H), x, c.start)
    Systems.evaluate_and_jacobian!(u, U, target(H), x, c.target)

    U .= (γ(H) * t) .* c.U .+ (1 - t) .* U
    u .= γ(H) .* c.u .- u

    nothing
end
