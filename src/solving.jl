module Solving

using Compat

import MultivariatePolynomials
const MP = MultivariatePolynomials

import ProgressMeter

import ..AffinePatches
import ..Endgame
import ..Homotopies
import ..PathTracking
import ..Parallel
import ..Problems
import ..PatchSwitching
import ..Systems

export Solver,
    solve

include("solving/path_result.jl")
include("solving/types.jl")
include("solving/solve.jl")

end
