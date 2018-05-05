module Homotopies

import ..Systems
import ..Systems: AbstractSystem, AbstractSystemCache
import Base: start

export AbstractHomotopy,
    StartTargetHomotopy,
    StraightLineHomotopy,
    γ,
    ParameterHomotopy,
    AbstractHomotopyCache,
    StartTargetHomotopyCache,
    nvariables,
    #start,
    target,
    cache,
    evaluate!, evaluate,
    jacobian!, jacobian,
    evaluate_and_jacobian!, evaluate_and_jacobian,
    dt!, dt,
    jacobian_and_dt!, jacobian_and_dt

export HomotopyWithCache

include("abstract_homotopy.jl")
include("straightline.jl")

end
