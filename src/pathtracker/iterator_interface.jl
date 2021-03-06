export track!, current_value

"""
    track!(pathtracker, x0, s_start, s_target)

Track a startvalue `x0` from `s_start` to `s_target` using the given `pathtracker`.

    track!(pathtracker)

Run the given `pathtracker`. You can use this in combination with [`setup_pathtracker!`](@ref).
"""
function track!(tracker::Pathtracker, startvalue::AbstractVector, start::Number=1.0, finish::Number=0.0, usehigh=false)
    setup_pathtracker!(tracker, startvalue, start, finish)
    if tracker.usehigh == false
        tracker.usehigh = usehigh
    end
    track!(tracker)
    tracker
end

function track!(tracker::Pathtracker)
    try
        start(tracker)
        is_done = done(tracker, 0)
        while !is_done
            next(tracker, 0)
            is_done = done(tracker, 0)
        end
    catch err
        if isa(err, Base.LinAlg.SingularException)
            tracker.status = :singular_exception
        else
            throw(err)
        end
    end
end


function Base.start(tracker::Pathtracker)
    precondition!(tracker, tracker.low, tracker.low.cache)

    if norm(evaluate(tracker.low.H, tracker.low.x, tracker.s, tracker.low.cfg)) > 0.1
        tracker.status = :invalid_startvalue
    end
    
    return 0
end

@inline function Base.next(tracker::Pathtracker, state)
    state += 1
    step!(tracker)

    tracker, state
end

@inline function Base.done(tracker::Pathtracker, state)
    if tracker.iter ≥ tracker.options.maxiters
        if tracker.options.verbose
            warn("Interrupted. Larger `maxiters` necessary.")
        end
        return true
    end
    if tracker.t ≤ 0.0
        return true
    end

    if tracker.status == :invalid_startvalue
        return true
    end

    false
end
Base.eltype(::T) where {T<:Pathtracker} = T


#TODO: Should a user need more hooks?
@inline function step!(tracker::Pathtracker)
    @inbounds pre_perform_step!(tracker)

    # this is implemented from the different algorithms / caches
    if tracker.usehigh
        @inbounds perform_step!(tracker, tracker.high, tracker.high.cache)
    else
        @inbounds perform_step!(tracker, tracker.low, tracker.low.cache)
    end

    @inbounds post_perform_step!(tracker)
end

@inline function pre_perform_step!(tracker::Pathtracker)
    tracker.ds = tracker.steplength * tracker.sdiff
    tracker.snext = tracker.s + tracker.ds
    tracker.iter += 1
    nothing
end

@inline function post_perform_step!(tracker::Pathtracker)
    if tracker.step_sucessfull
        tracker.t -= tracker.steplength

        tracker.consecutive_successfull_steps += 1
        if tracker.consecutive_successfull_steps ==
           tracker.options.consecutive_successfull_steps_until_steplength_increase
            tracker.steplength *= tracker.options.steplength_increase_factor
            tracker.consecutive_successfull_steps = 0
        end

        tracker.s = tracker.snext

        copy_xnext_to_x!(tracker)
    else
        tracker.consecutive_successfull_steps = 0
        tracker.steplength *= tracker.options.steplength_decrease_factor
    end
    tracker.steplength = min(tracker.steplength, tracker.t)
    nothing
end

@inline function copy_xnext_to_x!(tracker::Pathtracker)
    if tracker.usehigh
        copy!(tracker.high.x, tracker.high.xnext)
    else
        copy!(tracker.low.x, tracker.low.xnext)
    end
    nothing
end


is_projective_tracker(t::Pathtracker) = is_projective(t.alg)

"""
    current_value(pathtracker)

Get the current value of the pathtracker.
"""
function current_value(p::Pathtracker)
    if tracker.usehigh
        p.low.x .= p.high.x
    end
    copy(p.low.x)
end

# optional methods
precondition!(tracker, low, cache) = nothing

setup_workers(cache) = nothing
