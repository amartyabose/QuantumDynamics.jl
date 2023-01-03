module Utilities

using DifferentialEquations

function trapezoid(x, y; discrete::Bool=false)
    if discrete
        return sum(y)
    end
    sumvar = zero(y[1])
    for (a, b) in zip(y[2:end], y)
        sumvar += a + b
    end
    sumvar / 2 * (x[2] - x[1])
end

@inline function commutator(A, B)
    A * B - B * A
end

"""
    unhash_path(path_num::Int, ntimes::Int, sdim::Int)
Construct a path for a system with `sdim` dimensions, corresponding to the number `path_num`, with `ntimes` time steps.
"""
function unhash_path(path_num::Int, ntimes::Int, sdim::Int)
    path_num -= 1
    states = zeros(UInt8, ntimes + 1)
    for j in 1:ntimes+1
        @inbounds states[j] = path_num % sdim
        path_num = path_num ÷ sdim
    end
    states .+ 1
end

"""
    apply_propagator(; propagators, ρ0, ntimes, dt)
Apply a series of `ntimes` propagators to an initial reduced density matrix `ρ0` and return the result as a tuple of `(time, ρs)`.
"""
function apply_propagator(; propagators, ρ0, ntimes, dt)
    sdim = size(ρ0, 1)
    ρs = zeros(ComplexF64, ntimes + 1, sdim, sdim)
    @inbounds ρs[1, :, :] = ρ0
    ρvec = collect(Iterators.flatten(ρ0))
    for j = 1:ntimes
        @inbounds ρs[j+1, :, :] = reshape(propagators[j, :, :] * ρvec, (sdim, sdim))
    end
    0:dt:ntimes*dt, ρs
end

"""
ExternalField provides an abstract interface for encoding an external field, `V(t)`, interacting with the system through the operator, `coupling_op`.
"""
struct ExternalField
    V :: Function
    coupling_op :: Matrix{ComplexF64}
end

"""
Abstract type for encoding all the method specific numerical parameters.
"""
abstract type ExtraArgs end

"""
Extra parameters for solving differential equations. Currently has a threshold for magnitude-based filtering. The default values are:
    reltol = 1e-10
    abstol = 1e-10
    solver = Tsit5()
"""
struct DiffEqArgs <: Utilities.ExtraArgs
    reltol :: Float64
    abstol :: Float64
    solver
end
DiffEqArgs(; reltol=1e-10, abstol=1e-10, solver=Tsit5()) = DiffEqArgs(reltol, abstol, solver)

end