using LinearAlgebra

const statetostate_references = """
- A. Bose and P. L. Walters, Impact of Solvent on State-to-State Population Transport in Multistate Systems Using Coherences, J. Chem. Theory Comput. 19(15), 4828-4836 (2023).
- D. Sharma and A. Bose, Non-Hermitian State-to-State Analysis of Transport in Aggregates with Multiple Endpoints, J. Chem. Theory Comput. 21(12), 5858-5866 (2025).
- D. Sharma and A. Bose, Routes of Transport in the Path Integral Lindblad Dynamics through State-to-State Analysis, arXiv:2512.09362. https://doi.org/10.48550/arXiv.2512.09362"""
 
"""
    ddt_hamiltonian_flows(; ρs::AbstractArray{<:Complex,3}, H0::AbstractMatrix{<:Complex})
Returns a time-series of time derivative of state-to-state flows (\\dot P^H_{j \\leftarrow k}(t)) in each state (j) for a Hermitian or non-Hermitian system Hamiltonian H0.
"""
function ddt_hamiltonian_flows(; ρs::AbstractArray{<:Complex,3}, H0::AbstractMatrix{<:Complex})
    dim = size(H0, 1)
    N = size(ρs, 1)
    ddt_H0_flows = zeros(eltype(ρs), dim, N, dim)

    for j in 1:dim
        for n in 1:N
            for k in 1:dim
                ddt_H0_flows[j,n,k] = 1im * (ρs[n,j,k] * (H0')[k,j] - H0[j,k] * ρs[n,k,j])
            end
        end
    end

    ddt_H0_flows
end

kronecker_delta(i, j) = i == j ? 1 : 0

"""
    elementary_lindblad_states(l::Matrix{ComplexF64})
Returns the initial (\\ket{i_n}) and final states (\\ket{f_n}) affected by the action of a Lindblad jump operator l (=\\sum_n c_n\\dyad{f_n}{i_n}).
"""
function elementary_lindblad_states(l::Matrix{ComplexF64})
    ins = []
    fns = []
    dim = size(l, 1)
    for i in 1:dim
        for j in 1:dim
            if l[i,j] != (0.0 + 0.0im)
                @assert i in fns "Valid Lindblad jump operators restricted to forms where two different initial states do not map to the same final state. See reference: $(println("- D. Sharma and A. Bose, Routes of Transport in the Path Integral Lindblad Dynamics through State-to-State Analysis, arXiv:2512.09362. https://doi.org/10.48550/arXiv.2512.09362"))"
                push!(fns, i)
                push!(ins, j)
            end
        end
    end

    ins, fns
end

"""
    ddt_lindblad_flows(; ρs::AbstractArray{<:Complex,3}, L::Vector{Matrix{ComplexF64}})
Returns a time-series of time derivative of flows arising due to Lindblad pump(s) and/or drain(s) L for each state.
"""
function ddt_lindblad_flows(; ρs::AbstractArray{<:Complex,3}, L::Vector{Matrix{ComplexF64}})
    dim = size(ρs, 2)
    N = size(ρs, 1)
    ddt_L_flows = zeros(eltype(ρs), dim, N, dim)

    for j in 1:dim
        for n in 1:N
            for k in 1:dim
                for l in L
                    ins, fns = elementary_lindblad_states(l)
                    for i in eachindex(fns)
                        ddt_L_flows[j,n,k] += (l[fns[i],ins[i]]) ^ 2 * ρs[n,ins[i],ins[i]] * (kronecker_delta(j,fns[i]) * kronecker_delta(k, ins[i]) - kronecker_delta(j,ins[i]) * kronecker_delta(k, fns[i]))
                    end
                end
            end
        end
    end

    ddt_L_flows
end

"""
    statetostate(; t::AbstractArray{<:Real,1}, ρs::AbstractArray{<:Complex,3}, H0::AbstractMatrix{<:Complex}, L::Union{Nothing,Vector{Matrix{ComplexF64}}}=nothing)
Returns a time-series of total state-to-state flows (P_{j \\leftarrow k}(t)) (and its time derivative) in each state (j) which is the sum of flows arising due to Hamiltonian H0 and Lindbald jump operators L (if any).
"""
function statetostate(; t::AbstractArray{<:Real,1}, ρs::AbstractArray{<:Complex,3}, H0::AbstractMatrix{<:Complex}, L::Union{Nothing,Vector{Matrix{ComplexF64}}}=nothing)
    ddt_H0_flows = ddt_hamiltonian_flows(; ρs, H0)
    ddt_L_flows = zero(ddt_H0_flows)
    if !isnothing(L)
        ddt_L_flows = ddt_lindblad_flows(; ρs, L)
    end

    ddt_flows = ddt_H0_flows .+ ddt_L_flows

    flows = zero(ddt_flows)
    dim = size(H0, 1)
    N = length(t)

    for j in 1:dim
        for n in 1:N
            for k in 1:dim
                flows[j,n,k] = Utilities.trapezoid(t[1:n],ddt_flows[j,1:n,k])
            end
        end
    end

    ddt_flows, flows
end