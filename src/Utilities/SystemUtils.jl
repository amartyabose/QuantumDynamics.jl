using LinearAlgebra

"""
    create_nn_hamiltonian(; site_energies::AbstractVector{AbstractFloat}, couplings::AbstractVector{AbstractFloat}, periodic::Bool)
Creates a nearest neighbour Hamiltonian with the given `site_energies` and `couplings`. Periodic boundary conditions can also be used by passing `true` into the `periodic` argument.
"""
@inline function create_nn_hamiltonian(; site_energies::AbstractVector{<:AbstractFloat}, couplings::AbstractVector{<:AbstractFloat}, periodic::Bool)
    H = Array{Complex{real(eltype(site_energies))}}(diagm(0 => site_energies, 1 => couplings, -1 => couplings))
    if periodic
        H[1, end] += couplings
        H[end, 1] += couplings
    end
    H
end

"""
    create_tls_hamiltonian(; ϵ::AbstractFloat, Δ::AbstractFloat)

Creates a two-level system Hamiltonian:

``H = \\frac{ϵ}{2}σ_z - \\frac{Δ}{2}σ_x``

"""
create_tls_hamiltonian(; ϵ::AbstractFloat, Δ::AbstractFloat) = Array{Complex{typeof(ϵ)}}([ϵ/2 -Δ/2; -Δ/2 -ϵ/2])

"""
    calculate_Liouvillian(Hamiltonian::AbstractMatrix{Complex})
Returns the Liouvillian corresponding to the given Hamiltonian.
"""
function calculate_Liouvillian(Hamiltonian::AbstractMatrix{<:Complex})
    n = size(Hamiltonian, 1)
    identity_mat = Matrix{Complex{real(eltype(Hamiltonian))}}(I, n, n)
    kron(Hamiltonian, identity_mat) - kron(identity_mat, conj(Hamiltonian))
end

function calculate_Liouvillian(Hamiltonian::ITensor, sites)
    Hlat = swapinds(conj(Hamiltonian), [sites[1], sites[1]'], [sites[end], sites[end]'])
    Hamiltonian * delta(sites[end]', sites[end]) - delta(sites[1]', sites[1]) * Hlat
end

"""
ExternalField provides an abstract interface for encoding an external field, `V(t)`, interacting with the system through the operator, `coupling_op`.
"""
struct ExternalField
    V::Function
    coupling_op::Matrix{Complex}
end

"""
    apply_propagator(; propagators, ρ0, ntimes, dt)
Apply a series of `ntimes` propagators to an initial reduced density matrix `ρ0` and return the result as a tuple of `(time, ρs)`.
"""
function apply_propagator(; propagators, ρ0, ntimes, dt)
    sdim = size(ρ0, 1)
    ρs = zeros(eltype(propagators), ntimes + 1, sdim, sdim)
    @inbounds ρs[1, :, :] = ρ0
    ρvec = density_matrix_to_vector(ρ0)
    for j = 1:ntimes
        @inbounds ρs[j+1, :, :] = density_matrix_vector_to_matrix(propagators[j, :, :] * ρvec)
    end
    0:dt:ntimes*dt, ρs
end

"""
    commutator(A, B)
Returns the commutator A and B: AB - BA.
"""
function commutator(A, B)
    A * B .- B * A
end

"""
    density_matrix_to_vector(ρ::AbstractMatrix{<:Complex})
Returns the vector representation of the density matrix `ρ` compatible with the forward-backward propagators.
"""
density_matrix_to_vector(ρ::AbstractMatrix{<:Complex}) = collect(Iterators.flatten(transpose(ρ)))

"""
    density_matrix_vector_to_matrix(ρvec::AbstractVector{<:Complex})
Returns the matrix form of the vector `ρvec`.
"""
function density_matrix_vector_to_matrix(ρvec::AbstractVector{<:Complex})
    nsites = isqrt(length(ρvec))
    transpose(reshape(ρvec, (nsites, nsites)))
end

evaluate_observable(ρs::AbstractArray{ComplexF64, 3}, obs::AbstractMatrix{ComplexF64}) = [tr(ρs[j, :, :] * obs) for j in Axes(ρs, 1)]