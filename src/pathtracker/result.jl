export PathtrackerResult, solution

struct PathtrackerResult{T}
    returncode::Symbol
    solution::Vector{T}
    startvalue::Vector{T}
    residual::Float64
    iterations::Int
    angle_to_infinity::Float64

    # Extended analysis
    newton_residual::Float64
    condition_number::Float64
end

"""
    solution(pathtracker)

Get `(returncode, solution)` from `pathtracker`. This is more lightwheight than a
`PathtrackerResult`.
"""
@inline function solution(tracker::Pathtracker{Low}) where {Low}
    if tracker.iter ≥ tracker.options.maxiters
        returncode = :max_iterations
    elseif tracker.hit_singular_exception
        returncode = :singularity
    else
        returncode = :success
    end
    if tracker.usehigh
        sol = convert.(Complex{Low}, tracker.high.x)
    else
        sol = copy(tracker.low.x)
    end
    returncode, sol
end

function PathtrackerResult(tracker::Pathtracker{Low}, extended_analysis=true) where Low
    @unpack H, cfg = tracker.low
    returncode, sol = solution(tracker)

    res = evaluate(H, sol, tracker.s)
    residual = convert(Float64, norm(res))

    if extended_analysis
        jacobian = Homotopy.jacobian(H, sol, tracker.s, cfg)

        newton_residual = convert(Float64, norm(jacobian \ res))
        condition_number =  convert(Float64, Homotopy.κ(H, sol, 0.0, cfg))
    else
        newton_residual = NaN
        condition_number = NaN
    end

    if is_projective(tracker.alg)
        angle_to_infinity = convert(Float64, abs(first(sol)))
    else
        angle_to_infinity = 1.0
    end


    PathtrackerResult(
        returncode,
        sol,
        copy(tracker.startvalue),
        residual,
        tracker.iter,
        angle_to_infinity,
        newton_residual,
        condition_number)
end