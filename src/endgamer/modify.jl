export setup_endgamer!


"""
    setup_endgamer!(endgamer, x, R)

Setup `endgamer` to play the endgame starting from `x` at time `R`.
"""
function setup_endgamer!(endgamer::Endgamer, x, R)
    endgamer.R = R
    endgamer.iter = 0
    endgamer.windingnumber = 1
    endgamer.status = NotStarted
    endgamer.failurecode = :default_failure_code

    endgamer.startvalue .= x
    empty!(endgamer.xs)
    push!(endgamer.xs, x)

    empty!(endgamer.predictions)
end
