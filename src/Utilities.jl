module Utilities

using Combinatorics
using LinearAlgebra
using OrdinaryDiffEq
using ITensors

"""
    get_BLAS_implementation()
Reports the BLAS implementation under use. The default implementation used by Julia is OpenBlas. MKL is supported through the external package MKL.jl, which needs to be installed and loaded before the loading of QuantumDynamics.jl
"""
get_BLAS_implementation() = BLAS.get_config()

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
    states = zeros(Int, ntimes + 1)
    for j in 1:ntimes+1
        @inbounds states[j] = path_num % sdim
        path_num = path_num ÷ sdim
    end
    states .+ 1
end

function get_blip_starting_path(ntimes::Int, sdim::Int, nblips::Int, max::Int)
    if ntimes == 0
        return Vector{Vector{Int}}([])
    end
    if nblips == 0
        starting_paths = Vector{Vector{Int}}()
        push!(starting_paths, repeat([1], ntimes + 1))
        return starting_paths
    end
    starting_paths = Vector{Vector{Int}}()
    for l = 2:max
        if ntimes > 1
            rest = get_blip_starting_path(ntimes - 1, sdim, nblips - 1, l)
            for path in rest
                push!(starting_paths, vcat(path, l))
            end
        else
            if nblips == 1
                push!(starting_paths, [1, l])
            elseif nblips == 2
                for l2 = 2:l
                    push!(starting_paths, [l2, l])
                end
            end
        end
    end
    starting_paths
end

"""
    unhash_path_blips(ntimes::Int, sdim::Int, nblips::Int)
Construct all the paths for a system with `sdim` dimensions with `ntimes` time steps and `nblips` blips.
"""
function unhash_path_blips(ntimes::Int, sdim::Int, nblips::Int)
    starting_paths = get_blip_starting_path(ntimes, sdim, nblips, sdim)
    answers = Vector{Vector{Int}}()
    for p in starting_paths
        append!(answers, multiset_permutations(p, ntimes + 1) |> collect)
    end
    answers
end

"""
    apply_propagator(; propagators, ρ0, ntimes, dt)
Apply a series of `ntimes` propagators to an initial reduced density matrix `ρ0` and return the result as a tuple of `(time, ρs)`.
"""
function apply_propagator(; propagators, ρ0, ntimes, dt)
    sdim = size(ρ0, 1)
    ρs = zeros(ComplexF64, ntimes + 1, sdim, sdim)
    @inbounds ρs[1, :, :] = ρ0
    ρvec = collect(Iterators.flatten(transpose(ρ0)))
    for j = 1:ntimes
        @inbounds ρs[j+1, :, :] = transpose(reshape(propagators[j, :, :] * ρvec, (sdim, sdim)))
    end
    0:dt:ntimes*dt, ρs
end

""" 
    convert_ITensor_to_matrix(tens, sinit, sterm)
Converts an ITensor with two indices to a matrix. The index `sinit` is mapped to the column and `sterm` is mapped to the row.
"""
function convert_ITensor_to_matrix(tens, sinit, sterm)
    matrix = zeros(ComplexF64, dim(sterm), dim(sinit))
    for j = 1:dim(sterm)
        for k = 1:dim(sinit)
            matrix[j, k] = tens[sinit=>k, sterm=>j]
        end
    end
    matrix
end

"""
ExternalField provides an abstract interface for encoding an external field, `V(t)`, interacting with the system through the operator, `coupling_op`.
"""
struct ExternalField
    V::Function
    coupling_op::Matrix{ComplexF64}
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
    reltol::Float64
    abstol::Float64
    solver
end
DiffEqArgs(; reltol=1e-10, abstol=1e-10, solver=Tsit5()) = DiffEqArgs(reltol, abstol, solver)

"""
    create_nn_hamiltonian(; site_energies, couplings, periodic::Bool)
Creates a nearest neighbour Hamiltonian with the given `site_energies` and `couplings`. Periodic boundary conditions can also be used by passing `true` into the `periodic` argument.
"""
@inline function create_nn_hamiltonian(; site_energies, couplings, periodic::Bool)
    H = Array{ComplexF64}(diagm(0 => site_energies, 1 => couplings, -1 => couplings))
    if periodic
        H[1, end] += couplings
        H[end, 1] += couplings
    end
    H
end

"""
    create_tls_hamiltonian(; ϵ, Δ)

Creates a two-level system Hamiltonian:

``H = \\frac{ϵ}{2}σ_z - \\frac{Δ}{2}σ_x``

"""
create_tls_hamiltonian(; ϵ, Δ) = Array{ComplexF64}([ϵ/2+0.0im -Δ/2; -Δ/2 -ϵ/2])

end