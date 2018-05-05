module Correctors

using ..Homotopies
using ..Utilities

export AbstractCorrector,
    AbstractCorrectorCache,
    Result,
    cache,
    correct!,
    Newton, NewtonCache

include("abstract_corrector.jl")
include("newton.jl")

end
