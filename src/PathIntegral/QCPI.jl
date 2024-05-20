module QCPI

using ..EtaCoefficients, ..SpectralDensities, ..Solvents, ..Propagators, ..Utilities

"""
    propagate(; Hamiltonian::Matrix{ComplexF64}, Jw::SpectralDensities.SpectralDensity, solvent::Solvents.Solvent, ρ0::Matrix{ComplexF64}, classical_dt::Real, dt::Real, ntimes::Int, kmax::Int, reference_choice::String, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], verbose::Bool=false)

Use QCPI to propagate an initial density matrix, ρ0, under a given Hamiltonian with a solvent that is described by `solvent` and a corresponding spectral density `Jw`.
"""
function propagate(; Hamiltonian::Matrix{ComplexF64}, Jw::AbstractVector{<:SpectralDensities.SpectralDensity}, solvent::Solvents.Solvent, ρ0::Matrix{ComplexF64}, classical_dt::Real, dt::Real, ntimes::Int, kmax::Int, reference_choice::String, path_integral_routine, extraargs::Utilities.ExtraArgs, eacp=false, svec=[1.0 -1.0], verbose::Bool=false)
    sdim = size(Hamiltonian, 1)
    ρs = zeros(ComplexF64, ntimes + 1, sdim, sdim)
    time = 0:dt:ntimes*dt
    @assert length(Jw) == size(svec, 1)
    if kmax > ntimes
        kmax = ntimes + 2
    end
    ζ = [EtaCoefficients.calculate_ζ(jw; dt, kmax) for jw in Jw]
    for (i, ps) in enumerate(solvent)
        if verbose
            @info "Initial condition number $i started."
        end
        # fbU, srefs = Propagators.calculate_reference_propagators(; Hamiltonian, solvent, ps, ρ0, classical_dt, dt, ntimes, reference_choice)
        # if verbose
        #     @info "Propagators calculated."
        # end
        _, ρ = path_integral_routine(; β=solvent.β, Hamiltonian, solvent, ps, η=ζ, ρ0, dt, classical_dt, ntimes, kmax, extraargs, svec, reference_choice, verbose)
        ρs .+= ρ
        if verbose
            @info "Initial condition number $i done."
        end
    end
    time, ρs ./ solvent.num_samples
end

end
