module Propagators

using LinearAlgebra
using ..Solvents

function calculate_reference_propagators(; Hamiltonian::Matrix{ComplexF64}, solvent::Solvents.Solvent, classical_dt::Float64, dt::Float64, ntimes=1)
    nsys = size(Hamiltonian, 1)
    U = zeros(ComplexF64, ntimes, nsys^2, nsys^2)
    nclasstimes = trunc(Int, dt ÷ classical_dt)
    classical_dt = dt / nclasstimes
    ntotal_classical_steps = nclasstimes * ntimes
    for ps in solvent
        energy = Solvents.propagate_trajectory(solvent, ps, classical_dt, ntotal_classical_steps)
        last = 0
        for t = 1:ntimes
            Utmp = exp(-1im * (Hamiltonian .+ diagm(energy[t+last,:])) * classical_dt / 2.0)
            for j = 1:nclasstimes-1
                Utmp = exp(-1im * (Hamiltonian .+ diagm(energy[t+last+j,:])) * classical_dt) * Utmp
            end
            Utmp = exp(-1im * (Hamiltonian .+ diagm(energy[t+last+nclasstimes,:])) * classical_dt / 2.0) * Utmp
            last += nclasstimes-1
            U[t,:,:] .+= kron(Utmp, Utmp')
        end
    end
    U / solvent.num_samples
end

function calculate_bare_propagators(; Hamiltonian::Matrix{ComplexF64}, dt::Float64, ntimes=1)
    nsys = size(Hamiltonian, 1)
    U = zeros(ComplexF64, ntimes, nsys^2, nsys^2)
    for t = 1:ntimes
        U[t, :, :] .= kron(exp(-1im * Hamiltonian * dt), exp(1im * Hamiltonian * dt))
    end
    return U
end

end
