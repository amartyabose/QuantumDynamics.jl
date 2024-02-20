module Propagators

using LinearAlgebra
using Distributions
using ..Solvents, ..Utilities

make_fbpropagator(U) = kron(U, conj(U))
make_fbpropagator(U, Udag) = kron(U, transpose(Udag))

function calculate_reference_propagators(; Hamiltonian::AbstractMatrix{<:Complex}, solvent::Solvents.Solvent, ps::Solvents.PhaseSpace, classical_dt::AbstractFloat, dt::AbstractFloat, ref_pos=nothing, ntimes=1)
    nsys = size(Hamiltonian, 1)
    U = zeros(eltype(Hamiltonian), ntimes, nsys^2, nsys^2)
    nclasstimes = trunc(Int, dt ÷ classical_dt)
    classical_dt = dt / nclasstimes
    ntotal_classical_steps = nclasstimes * ntimes
    ref_pos_mod = zeros(ntotal_classical_steps)
    if !isnothing(ref_pos)
        count = 1
        for j = 1:ntimes
            for k = 1:nclasstimes
                @inbounds ref_pos_mod[count] = ref_pos[j]
                count += 1
            end
        end
    end
    energy, _ = Solvents.propagate_trajectory(solvent, ps, classical_dt, ntotal_classical_steps, ref_pos_mod)
    Href = copy(Hamiltonian)
    last = 0
    diaginds = diagind(Href)
    eye = Matrix{eltype(Hamiltonian)}(I, nsys, nsys)
    for t = 1:ntimes
        @inbounds Href .= Hamiltonian
        @inbounds Href[diaginds] .+= energy[t+last, :]
        @inbounds Href .-= eye * tr(Href) / nsys
        Utmp = exp(-1im * Href * classical_dt / 2.0)
        for j = 1:nclasstimes-1
            @inbounds Href .= Hamiltonian
            @inbounds Href[diaginds] .+= energy[t+last+j, :]
            @inbounds Href .-= eye * tr(Href) / nsys
            Utmp = exp(-1im * Href * classical_dt) * Utmp
        end
        @inbounds Href .= Hamiltonian
        @inbounds Href[diaginds] .+= energy[t+last+nclasstimes, :]
        @inbounds Href .-= eye * tr(Href) / nsys
        Utmp = exp(-1im * Href * classical_dt / 2.0) * Utmp
        last += nclasstimes - 1
        @inbounds U[t, :, :] .+= make_fbpropagator(Utmp)
    end
    U
end

function calculate_average_reference_propagators(; Hamiltonian::AbstractMatrix{<:Complex}, solvent::Solvents.Solvent, classical_dt::AbstractFloat, dt::AbstractFloat, ntimes=1, verbose=false)
    nsys = size(Hamiltonian, 1)
    elem_type = eltype(Hamiltonian)
    U = zeros(elem_type, ntimes, nsys^2, nsys^2)
    Ucum = zeros(elem_type, ntimes, nsys^2, nsys^2)
    Utmp = zeros(elem_type, ntimes, nsys^2, nsys^2)
    ref_pos = zeros(ntimes + 1)
    for (i, ps) in enumerate(solvent)
        @inbounds Utmp .= calculate_reference_propagators(; Hamiltonian, solvent, ps, classical_dt, dt, ref_pos, ntimes)
        @inbounds Ucum[1, :, :] .= Utmp[1, :, :]
        for j = 2:ntimes
            @inbounds Ucum[j, :, :] .= Ucum[j-1, :, :] * Utmp[j, :, :]
        end
        @inbounds U .+= Ucum
        if verbose
            @info "Phase space point $(i) of $(solvent.num_samples)"
        end
    end
    U ./ solvent.num_samples
end

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
