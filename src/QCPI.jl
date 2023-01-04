module QCPI

using ..SpectralDensities, ..Solvents, ..Propagators, ..Utilities

"""
    propagate(; Hamiltonian::Matrix{ComplexF64}, Jw::SpectralDensities.SpectralDensity, solvent::Solvents.Solvent, ρ0::Matrix{ComplexF64}, classical_dt::Real, dt::Real, ntimes::Int, kmax::Int, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], verbose::Bool=false)

Use QCPI to propagate an initial density matrix, ρ0, under a given Hamiltonian with a solvent that is described by `solvent` and a corresponding spectral density `Jw`.
"""
function propagate(; Hamiltonian::Matrix{ComplexF64}, Jw::SpectralDensities.SpectralDensity, solvent::Solvents.Solvent, ρ0::Matrix{ComplexF64}, classical_dt::Real, dt::Real, ntimes::Int, kmax::Int, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], verbose::Bool=false)
    sdim = size(Hamiltonian, 1)
    ρs = zeros(ComplexF64, ntimes+1, sdim, sdim)
    time = 0:dt:ntimes*dt
    for ps in solvent
        fbU = Propagators.calculate_reference_propagators(; Hamiltonian, solvent, ps, classical_dt, dt, ntimes)
        _, ρ = path_integral_routine(; fbU, Jw=[Jw], β=solvent.β, ρ0, dt, ntimes, kmax, extraargs, svec, reference_prop=true, verbose)
        ρs .+= ρ
    end
    time, ρs ./ solvent.num_samples
end

end