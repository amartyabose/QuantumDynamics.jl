module Propagators

using LinearAlgebra
using Distributions
using FLoops

using ITensors, ITensorMPS

using ..Solvents, ..Utilities

make_fbpropagator(U) = kron(U, conj(U))
make_fbpropagator(U, Udag) = kron(U, transpose(Udag))

"""
    ReferenceChoice{choice}
This empty struct holds the particular choice of algorithm for selecting reference forces. `choice` currently can be one of `fixed` or `ehrenfest`.
"""
struct ReferenceChoice{choice} end
ReferenceChoice(c) = ReferenceChoice{Symbol(c)}
macro ReferenceChoice_str(c)
    return :(ReferenceChoice{$(Expr(:quote, Symbol(c)))})
end

get_reference(::ReferenceChoice"fixed", ρ, solvent::Solvents.Solvent) = zeros(solvent.num_baths)
get_reference(::ReferenceChoice"ehrenfest", ρ, solvent::Solvents.Solvent) = [real(tr(ρ * diagm(solvent.sys_op[j, :]))) for j = 1:solvent.num_baths]
get_reference(::ReferenceChoice"max_prob", ρ, solvent::Solvents.Solvent) = [solvent.sys_op[j, findmax(real.(diag(ρ)))[2]] for j = 1:solvent.num_baths]
function get_reference(::ReferenceChoice"dcsh", ρ, solvent::Solvents.Solvent)
    vecd = real.(diag(ρ))
    pos = findfirst(x -> x>1.0, vecd)
    if !isnothing(pos)
        [solvent.sys_op[j, pos] for j = 1:solvent.num_baths]
    else
        vecd ./= sum(vecd)
        distrib = Distributions.Categorical(vecd)
        r = rand(distrib)
        [solvent.sys_op[j, r] for j = 1:solvent.num_baths]
    end
end
function get_reference(::ReferenceChoice"ehrenfest_surface", ρ, solvent::Solvents.Solvent)
    vals = [real(tr(ρ * diagm(solvent.sys_op[j, :]))) for j = 1:solvent.num_baths]
    ref = zeros(length(vals))
    for j in eachindex(vals)
        ref[j] = solvent.sys_op[j, findmin(abs.(solvent.sys_op[j, :] .- vals[j]))[2]]
    end
    ref
end

function calculate_reference_propagators(; Hamiltonian::AbstractMatrix{<:Complex}, solvent::Solvents.Solvent, ps::Solvents.PhaseSpace, classical_dt::AbstractFloat, dt::AbstractFloat, ρ0::AbstractMatrix{<:Complex}, ntimes=1, reference_choice::String="ehrenfest")
    nsys = size(Hamiltonian, 1)
    U = zeros(eltype(Hamiltonian), ntimes, nsys^2, nsys^2)
    nclasstimes = trunc(Int, dt ÷ classical_dt)
    classical_dt = dt / nclasstimes
    Href = copy(Hamiltonian)
    diaginds = diagind(Href)
    eye = Matrix{eltype(Hamiltonian)}(I, nsys, nsys)
    ρ = copy(ρ0)
    srefs = zeros(ntimes, length(solvent.num_baths))
    for t = 1:ntimes
        reference = get_reference(ReferenceChoice(reference_choice)(), ρ, solvent)#, Hamiltonian, dt, ps)
        energy, ps = Solvents.propagate_trajectory(solvent, ps, classical_dt, nclasstimes, reference)
        @inbounds srefs[t, :] .= reference
        @inbounds Href .= Hamiltonian
        @inbounds Href[diaginds] .+= energy[1, :]
        @inbounds Href .-= eye * tr(Href) / nsys
        Utmp = exp(-1im * Href * classical_dt / 2.0)
        for j = 1:nclasstimes-1
            @inbounds Href .= Hamiltonian
            @inbounds Href[diaginds] .+= energy[1+j, :]
            @inbounds Href .-= eye * tr(Href) / nsys
            Utmp = exp(-1im * Href * classical_dt) * Utmp
        end
        @inbounds Href .= Hamiltonian
        @inbounds Href[diaginds] .+= energy[1, :]
        @inbounds Href .-= eye * tr(Href) / nsys
        Utmp = exp(-1im * Href * classical_dt / 2.0) * Utmp
        ρ = Utmp * ρ * Utmp'
        @inbounds U[t, :, :] .+= make_fbpropagator(Utmp)
    end
    U, srefs, ps
end

function calculate_average_reference_propagators(; Hamiltonian::AbstractMatrix{<:Complex}, solvent::Solvents.Solvent, classical_dt::AbstractFloat, dt::AbstractFloat, ρ0, ntimes=1, verbose=false, reference_choice::String="ehrenfest")
    verbose && @info "Using $(reference_choice) algorithm for choosing the system state for the reference trajectory."
    nsys = size(Hamiltonian, 1)
    elem_type = eltype(Hamiltonian)
    U = zeros(elem_type, ntimes, nsys^2, nsys^2)
    for (i, ps) in enumerate(solvent)
        Utmp, _, _ = calculate_reference_propagators(; Hamiltonian, solvent, ps, classical_dt, dt, ρ0, ntimes, reference_choice)
        for j = 2:size(Utmp, 1)
            @inbounds Utmp[j, :, :] = Utmp[j-1, :, :] * Utmp[j, :, :]
        end
        @inbounds U .+= Utmp
        if verbose# && (i==1 || trunc(Int64, i / length(solvent) * 100) ÷ 10 > trunc(Int64, (i-1) / length(solvent) * 100) ÷ 10)
            @info "Phase space point $(i) of $(solvent.num_samples)"
        end
    end
    U ./ solvent.num_samples
end

function calculate_average_reference_propagators_parallel(; Hamiltonian::AbstractMatrix{<:Complex}, solvent::Solvents.Solvent, classical_dt::AbstractFloat, dt::AbstractFloat, ρ0, ntimes=1, verbose=false, reference_choice::String="ehrenfest")
    @info "Using $(reference_choice) algorithm for choosing the system state for the reference trajectory."
    nsys = size(Hamiltonian, 1)
    elem_type = eltype(Hamiltonian)
    U = zeros(elem_type, ntimes, nsys^2, nsys^2)
    pspoints = collect(solvent)
    @floop for ps in pspoints
        Utmp, _, _ = calculate_reference_propagators(; Hamiltonian, solvent, ps, classical_dt, dt, ρ0, ntimes, reference_choice)
        for j = 2:size(Utmp, 1)
            @inbounds Utmp[j, :, :] = Utmp[j-1, :, :] * Utmp[j, :, :]
        end
        @reduce U .= zeros(elem_type, ntimes, nsys^2, nsys^2) .+ Utmp
    end
    U ./ solvent.num_samples
end

function calculate_average_reference_propagators_mps(; Hamiltonian::AbstractMatrix{<:Complex}, solvent::Solvents.Solvent, classical_dt::AbstractFloat, dt::AbstractFloat, ref_pos::Union{Nothing,Vector{Float64}}=nothing, ntimes=1, verbose=false, reference_choice::String="ehrenfest", extraargs::Utilities.TensorNetworkArgs=Utilities.TensorNetworkArgs(; algorithm="densitymatrix"))
    nsys2 = size(Hamiltonian, 1)^2
    sites = siteinds(nsys2, ntimes + 1)
    Us = Vector{MPS}([])
    for (i, ps) in enumerate(solvent)
        Utmp, _, _ = calculate_reference_propagators(; Hamiltonian, solvent, ps, classical_dt, dt, ref_pos, ntimes, reference_choice)
        Umps = [Utilities.build_path_amplitude_mps(Utmp[1, :, :], sites[1:2])]
        for j = 2:ntimes
            push!(Umps, Utilities.extend_path_amplitude_mps(Umps[end], Utmp[j, :, :], sites[j:j+1]))
        end
        if i == 1
            Us = Umps
        else
            for j = 1:ntimes
                Us[j] = add(Us[j], Umps[j]; cutoff=extraargs.cutoff, maxdim=extraargs.maxdim, alg=extraargs.algorithm)
            end
        end
        if verbose
            @info "Phase space point $(i) of $(solvent.num_samples). Max bond dim = $(maximum.(linkdims.(Us)))"
        end
    end
    Us ./ solvent.num_samples
end

"""
    calculate_bare_propagators(; Hamiltonian::AbstractMatrix{<:Complex}, dt::AbstractFloat, ntimes=1, external_fields::Union{Nothing,Vector{Utilities.ExternalField}}=nothing, forward_backward=true, L::Union{Nothing, Vector{Matrix{ComplexF64}}}=nothing)
This function calculates the bare propagators for a given system.

Arguments:
- `Hamiltonian`: the system is fundamentally described by the Hamiltonian
- `external_fields` [Default: nothing]: an optional listing of external fields that are acting on the system
- `forward_backward` [Default: true]: do you want a forward-backward propagator or just the forward propagator?
- `L` [Default: nothing]: an optional listing of Lindblad jump operators acting on the system (Acts only when `forward_backward = true`)
- `dt`: time-step corresponding to the bare propagator calculation
- `ntimes` [Default: 1]: number of time-steps for which the bare propagator should be calculated

Returns the `ntimes` propagators.
"""
function calculate_bare_propagators(; Hamiltonian::AbstractMatrix{<:Complex}, dt::AbstractFloat, ntimes=1, external_fields::Union{Nothing,Vector{Utilities.ExternalField}}=nothing, forward_backward=true, L::Union{Nothing, Vector{Matrix{ComplexF64}}}=nothing)
    nsys = size(Hamiltonian, 1)
    U = forward_backward ? zeros(eltype(Hamiltonian), ntimes, nsys^2, nsys^2) : zeros(eltype(Hamiltonian), ntimes, nsys, nsys)
    if isnothing(external_fields)
        if forward_backward
            liouvillian = Utilities.calculate_Liouvillian(Hamiltonian)
            identity_mat = Matrix{Complex{real(eltype(Hamiltonian))}}(I, nsys, nsys)
            if !isnothing(L)
                for l in L
                    ldagl = l' * l
                    liouvillian += kron(l, conj(l)) - 0.5 * (kron(ldagl, identity_mat) + kron(identity_mat, transpose(ldagl)))
                end
            end
            Utmp = exp(liouvillian * dt)
            for t = 1:ntimes
                @inbounds U[t, :, :] .+= Utmp
            end
        else
            Utmp = exp(-1im * Hamiltonian * dt)
            for t = 1:ntimes
                @inbounds U[t, :, :] .+= Utmp
            end
        end
    else
        ndivs = 10000
        delt = dt / ndivs
        if forward_backward
            identity_mat = Matrix{Complex{real(eltype(Hamiltonian))}}(I, nsys, nsys)
            liouvillian_jump = zeros(nsys^2, nsys^2)
            if !isnothing(L)
                for l in L
                    ldagl = l' * l
                    liouvillian_jump += kron(l, conj(l)) - 0.5 * (kron(ldagl, identity_mat) + kron(identity_mat, conj(ldagl)))
                end
            end
            for t = 1:ntimes
                Utmp = Matrix{eltype(Hamiltonian)}(I, nsys^2, nsys^2)
                for j = 1:ndivs
                    H = copy(Hamiltonian)
                    for ef in external_fields
                        @inbounds H .+= ef.V(((t - 1) * ndivs + (j - 1)) * delt) .* ef.coupling_op
                    end
                    Utmp = exp((Utilities.calculate_Liouvillian(H) + liouvillian_jump) * delt) * Utmp
                end
                @inbounds U[t, :, :] .+= Utmp
            end
        else
            for t = 1:ntimes
                Utmp = Matrix{eltype(Hamiltonian)}(I, nsys, nsys)
                for j = 1:ndivs
                    H = copy(Hamiltonian)
                    for ef in external_fields
                        @inbounds H .+= ef.V(((t - 1) * ndivs + (j - 1)) * delt) .* ef.coupling_op
                    end
                    Utmp = exp(-1im * H * delt) * Utmp
                end
                @inbounds U[t, :, :] .+= Utmp
            end
        end
    end
    U
end

end
