module QCPI

using ProgressMeter
using ..SpectralDensities, ..Solvents, ..Propagators

function propagate(; Hamiltonian::Matrix{ComplexF64}, Jw::SpectralDensities.SpectralDensity, solvent::Solvents.Solvent, ρ0::Matrix{ComplexF64}, classical_dt::Real, dt::Real, ntimes::Int, kmax::Int, path_integral_routine, cutoff=0.0, svec=[1.0 -1.0], verbose::Bool=false)
    sdim = size(Hamiltonian, 1)
    ρs = zeros(ComplexF64, ntimes+1, sdim, sdim)
    time = 0:dt:ntimes*dt
    for ps in solvent
        fbU = Propagators.calculate_reference_propagators(; Hamiltonian, solvent, ps, classical_dt, dt, ntimes)
        _, ρ = path_integral_routine(; fbU, Jw=[Jw], β=solvent.β, ρ0, dt, ntimes, kmax, cutoff, svec, reference_prop=true, verbose)
        ρs .+= ρ
    end
    time, ρs ./ solvent.num_samples
end

end