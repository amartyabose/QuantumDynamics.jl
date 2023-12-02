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

"""
    trapezoid(x, y; discrete::Bool=false)
Returns the trapezoidal integration of y with respect to x. If discrete is set to `true`, then returns sum of y.
"""
trapezoid(x, y; discrete::Bool=false) = discrete ? sum(y) : sum((y[1:end-1] .+ y[2:end]) .* (x[2:end] .- x[1:end-1])) / 2

function fourier_transform(time, corr; full=true, unitary=false)
    dt = time[2] - time[1]
    ωmax = π / dt
    dω = π / time[end]
    ω = -ωmax:dω:ωmax |> collect
    spect = zeros(ComplexF64, length(ω))
    if full
        for (l, w) in enumerate(ω)
            spect[l] = Utilities.trapezoid(time, corr .* exp.(-1im * w * time) + conj.(corr) .* exp.(1im * w * time))
        end
    else
        for (l, w) in enumerate(ω)
            spect[l] = Utilities.trapezoid(time, corr .* exp.(-1im * w * time))
        end
    end

    spect ./= unitary ? sqrt(2π) : 1

    ω, spect
end

function inverse_fourier_transform(ω, spect; unitary=false)
    dω = ω[2] - ω[1]
    dt = π / ω[end]
    tmax = π / dω
    time = 0:dt:tmax
    data = zeros(ComplexF64, length(time))
    for (l, t) in enumerate(time)
        data[l] = Utilities.trapezoid(ω, spect .* exp.(1im * ω * t))
    end

    data ./= unitary ? sqrt(2π) : 2π
    time, data
end

"""
    commutator(A, B)
Returns the commutator A and B: AB - BA.
"""
function commutator(A, B)
    A * B .- B * A
end

"""
    unhash_path(path_num::Int, ntimes::Int, sdim::Int)
Construct a path for a system with `sdim` dimensions, corresponding to the number `path_num`, with `ntimes` time steps.
"""
function unhash_path(path_num::Int, ntimes::Int, sdim)
    path_num -= 1
    states = zeros(UInt8, ntimes + 1)
    for j in 1:ntimes+1
        @inbounds states[j] = path_num % sdim
        path_num = path_num ÷ sdim
    end
    states .+ 1
end

"""
    hash_path(states, sdim)
Returns the hashed location of a path for a system with `sdim` dimensions.
"""
function hash_path(states::Vector{UInt8}, sdim)
    factor = 1
    number = 0
    for s in states
        number += (s - 1) * factor
        factor *= sdim
    end
    number + 1
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
    calculate_Liouvillian(Hamiltonian::AbstractMatrix{ComplexF64})
Returns the Liouvillian corresponding to the given Hamiltonian.
"""
function calculate_Liouvillian(Hamiltonian::AbstractMatrix{ComplexF64})
    n = size(Hamiltonian, 1)
    identity_mat = Matrix{ComplexF64}(I, n, n)
    kron(Hamiltonian, identity_mat) - kron(identity_mat, conj(Hamiltonian))
end

"""
    calculate_Liouvillian(H::OpSum, sites)
Returns the forward-backward space combiner and the Liouvillian MPO corresponding to the Hamiltonian provided as an ITensor `OpSum`, to be built on the `sites`.
"""
function calculate_Liouvillian(H::OpSum, sites)
    dupsites = addtags(sites, "*")
    fbcombiner = forward_backward_combiner(sites, dupsites)

    Hamiltonian = MPO(H, sites)
    Hamiltonian_tilde = MPO(H, dupsites)
    id_MPO = identity_MPO(sites)
    id_MPO_tilde = identity_MPO(dupsites)
    liouvillian = Hamiltonian * id_MPO_tilde - id_MPO * Hamiltonian_tilde
    for j = 1:length(sites)
        liouvillian[j] = liouvillian[j] * fbcombiner[j] * fbcombiner[j]'
    end

    fbcombiner, liouvillian
end

forward_backward_combiner(sites, sites1=sites) = [combiner(s, s1; tags="FBSite") for (s, s1) in zip(sites, sites1)]

"""
    identity_MPO(sites)
Returns the identity MPO based on the given `sites`.
"""
function identity_MPO(sites)
    idMPO = MPO(sites)
    linds = linkinds(idMPO)
    sitedim = dims(sites)[1]
    for s = 1:sitedim
        idMPO[1][linds[1]=>1, sites[1]=>s, sites[1]'=>s] = 1
    end
    for j = 2:length(sites)-1
        sitedim = dims(sites)[j]
        for s = 1:sitedim
            idMPO[j][linds[j]=>1, linds[j-1]=>1, sites[j]=>s, sites[j]'=>s] = 1
        end
    end
    sitedim = dims(sites)[end]
    for s = 1:sitedim
        idMPO[end][linds[end]=>1, sites[end]=>s, sites[end]'=>s] = 1
    end
    idMPO
end

"""
    MPO_to_MPS(ρ::MPO, fbcombiner)
Convert a given MPO to an MPS by combining the two site indices on every tensor into a single one using the vector of combiner tensors.
"""
function MPO_to_MPS(ρ::MPO, fbcombiner)
    ρfb = deepcopy(ρ)
    for j = 1:length(ρ)
        ρfb[j] *= fbcombiner[j]
    end
    convert(MPS, ρfb)
end

"""
    MPO_to_MPS(ρ::MPO, fbcombiner)
Split a given MPS to an MPO by using the vector of combiner tensors.
"""
function MPS_to_MPO(ρ::MPS, fbcombiner)
    ρfb = deepcopy(ρ)
    for j = 1:length(ρ)
        ρfb[j] *= fbcombiner[j]
    end
    convert(MPO, ρfb)
end

"""
    ITensors.expect(ρ::MPO, ops; kwargs...)
Extends ITensors' `expect` function to handle density matrices in the form of MPOs.
"""
function ITensors.expect(ρ::MPO, ops::Tuple; kwargs...)
    ρtmp = deepcopy(ρ)
    N = length(ρ)
    s = Vector{Index{Int64}}()
    sstar = Vector{Index{Int64}}()
    for ρind in siteinds(ρ)
        for j in ρind
            if hastags(j, "*")
                push!(sstar, j)
            else
                push!(s, j)
            end
        end
    end
    # s = [(map(filter(x->!hastags(x,"*")), siteinds(ρ))...)...]
    # sstar = [(map(filter(x->hastags(x,"*")), siteinds(ρ))...)...]
    for j = 1:N
        swapinds!(ρtmp[j], sstar[j], s[j]')
    end

    if haskey(kwargs, :site_range)
        @warn "The `site_range` keyword arg. to `expect` is deprecated: use the keyword `sites` instead"
        sites = kwargs[:site_range]
    else
        sites = get(kwargs, :sites, 1:N)
    end

    site_range = (sites isa AbstractRange) ? sites : collect(sites)
    Ns = length(site_range)
    start_site = first(site_range)

    el_types = map(o -> ishermitian(op(o, s[start_site])) ? Float64 : ComplexF64, ops)

    ex = map((o, el_t) -> zeros(el_t, Ns), ops, el_types)
    for (entry, j) in enumerate(site_range)
        for (n, opname) in enumerate(ops)
            ans = j == 1 ? swapprime(ρtmp[j] * op(opname, s[j])', 2 => 1) * delta(s[j], s[j]') : ρtmp[1] * delta(s[1], s[1]')
            for k = 2:j-1
                ans *= ρtmp[k] * delta(s[k], s[k]')
            end
            if j != 1
                ans *= swapprime(ρtmp[j] * op(opname, s[j])', 2 => 1) * delta(s[j], s[j]')
            end
            for k = j+1:N
                ans *= ρtmp[k] * delta(s[k], s[k]')
            end
            val = scalar(ans)
            ex[n][entry] = (el_types[n] <: Real) ? real(val) : val
        end
    end

    if sites isa Number
        return map(arr -> arr[1], ex)
    end
    return ex
end

function ITensors.expect(psi::MPO, op::AbstractString; kwargs...)
    return first(expect(psi, (op,); kwargs...))
end

function ITensors.expect(psi::MPO, op::Matrix{<:Number}; kwargs...)
    return first(expect(psi, (op,); kwargs...))
end

function ITensors.expect(psi::MPO, op1::AbstractString, ops::AbstractString...; kwargs...)
    return expect(psi, (op1, ops...); kwargs...)
end

function ITensors.expect(psi::MPO, op1::Matrix{<:Number}, ops::Matrix{<:Number}...; kwargs...)
    return expect(psi, (op1, ops...); kwargs...)
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