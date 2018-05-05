__precompile__()

module HomotopyContinuation

    import DynamicPolynomials: @polyvar
    export @polyvar

    export AffinePatches,
        Correctors,
        Endgame,
        Homotopies,
        PathTracking,
        Predictors,
        Problems,
        ProjectiveVectors,
        Solving,
        StepLength,
        Systems,
        Utilities

    include("utilities.jl")
    include("parallel.jl")
    include("projective_vectors.jl")
    include("systems.jl")
    include("homotopies/_Homotopies.jl")
    include("problems.jl")
    include("predictors/_Predictors.jl")
    include("correctors/_Correctors.jl")
    include("prediction_correction.jl")
    include("affine_patches.jl")
    include("step_length.jl")

    include("path_tracking.jl")
    include("endgame.jl")

    include("solving.jl")
    include("solve.jl")
end #
