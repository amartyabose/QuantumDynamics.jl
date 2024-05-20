module Propagators

using LinearAlgebra
using Distributions
using FLoops

using ITensors

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

get_reference(::ReferenceChoice"fixed", ρ, solvent::Solvents.Solvent) = zeros(solvent.num_modes)
get_reference(::ReferenceChoice"ehrenfest", ρ, solvent::Solvents.Solvent) = [real(tr(ρ * diagm(solvent.sys_op[j, :]))) for j = 1:solvent.num_modes]
get_reference(::ReferenceChoice"max_prob", ρ, solvent::Solvents.Solvent) = [solvent.sys_op[j, findmax(real.(diag(ρ)))[2]] for j = 1:solvent.num_modes]
function get_reference(::ReferenceChoice"dcsh", ρ, solvent::Solvents.Solvent)
    vec = real.(diag(ρ))
    clamp!(vec, 0.0, 1.0)
    vec ./= sum(vec)
    distrib = Distributions.Categorical(vec)
    r = rand(distrib)
    [solvent.sys_op[j, r] for j = 1:solvent.num_modes]
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
    srefs = zeros(ntimes, length(solvent.init_points))
    for t = 1:ntimes
        reference = get_reference(ReferenceChoice(reference_choice)(), ρ, solvent)#, Hamiltonian, dt, ps)
        energy, _, ps = Solvents.propagate_trajectory(solvent, ps, classical_dt, nclasstimes, reference)
        @inbounds srefs[t, :] .= reference[solvent.init_points]
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
    nsys = size(Hamiltonian, 1)
    elem_type = eltype(Hamiltonian)
    U = zeros(elem_type, ntimes, nsys^2, nsys^2)
    for (i, ps) in enumerate(solvent)
        Utmp, _, _ = calculate_reference_propagators(; Hamiltonian, solvent, ps, classical_dt, dt, ρ0, ntimes, reference_choice)
        for j = 2:size(Utmp, 1)
            @inbounds Utmp[j, :, :] = Utmp[j-1, :, :] * Utmp[j, :, :]
        end
        @inbounds U .+= Utmp
        if verbose
            @info "Phase space point $(i) of $(solvent.num_samples)"
        end
    end
    U ./ solvent.num_samples
end

function calculate_average_reference_propagators_parallel(; Hamiltonian::AbstractMatrix{<:Complex}, solvent::Solvents.Solvent, classical_dt::AbstractFloat, dt::AbstractFloat, ρ0, ntimes=1, verbose=false, reference_choice::String="ehrenfest")
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
    calculate_bare_propagators(; Hamiltonian::AbstractMatrix{<:Complex}, dt::AbstractFloat, ntimes=1, external_fields::Union{Nothing,Vector{Utilities.ExternalField}}=nothing)
This function calculates the bare propagators for a given `Hamiltonian` and under the influence of the `external_fields` with a time-step of `dt` for `ntimes` time steps.
"""
function calculate_bare_propagators(; Hamiltonian::AbstractMatrix{<:Complex}, dt::AbstractFloat, ntimes=1, external_fields::Union{Nothing,Vector{Utilities.ExternalField}}=nothing)
    nsys = size(Hamiltonian, 1)
    U = zeros(eltype(Hamiltonian), ntimes, nsys^2, nsys^2)
    if isnothing(external_fields)
        Utmp = exp(-1im * Hamiltonian * dt)
        Utmpdag = exp(1im * Hamiltonian' * dt)
        for t = 1:ntimes
            @inbounds U[t, :, :] .+= make_fbpropagator(Utmp, Utmpdag)
        end
    else
        ndivs = 10000
        delt = dt / ndivs
        for t = 1:ntimes
            Utmp = Matrix{eltype(Hamiltonian)}(I, nsys, nsys)
            Utmpdag = Matrix{eltype(Hamiltonian)}(I, nsys, nsys)
            for j = 1:ndivs
                H = copy(Hamiltonian)
                for ef in external_fields
                    @inbounds H .+= ef.V(((t - 1) * ndivs + (j - 1)) * delt) .* ef.coupling_op
                end
                Utmp = exp(-1im * H * delt) * Utmp
                Utmpdag = exp(1im * H' * delt) * Utmpdag
            end
            @inbounds U[t, :, :] .+= make_fbpropagator(Utmp, Utmpdag)
        end
    end
    U
end

end
