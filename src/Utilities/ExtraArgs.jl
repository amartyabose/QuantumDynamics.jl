using OrdinaryDiffEq

"""
Abstract type for encoding all the method specific numerical parameters.
"""
abstract type ExtraArgs end

"""
Extra parameters for solving differential equations. Currently has a threshold for magnitude-based filtering. The default values are:
- reltol = 1e-10
- abstol = 1e-10
- solver = Tsit5()
"""
struct DiffEqArgs <: Utilities.ExtraArgs
    reltol::Float64
    abstol::Float64
    solver
end
DiffEqArgs(; reltol=1e-10, abstol=1e-10, solver=Tsit5()) = DiffEqArgs(reltol, abstol, solver)

"""
Extra parameters for tensor network algorithms. Currently has an SVD `cutoff`, a maximum bond dimension `maxdim`, and the contraction `algorithm`. The default values are as follows:
- cutoff = 1e-8
- maxdim = 500
- algorithm = "naive"

Other options for algorithm are "densitymatrix" and "fit".
"""
struct TensorNetworkArgs <: Utilities.ExtraArgs
    cutoff::AbstractFloat
    maxdim::Integer
    algorithm::String
end
TensorNetworkArgs(; cutoff=1e-8, maxdim=500, algorithm="naive") = TensorNetworkArgs(cutoff, maxdim, algorithm)
