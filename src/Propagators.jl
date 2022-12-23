module Propagators

using LinearAlgebra
using DifferentialEquations
using ..Solvents

function calculate_reference_propagators(; Hamiltonian::Matrix{ComplexF64}, solvent::Solvents.Solvent, ps::Solvents.PhaseSpace, classical_dt::Float64, dt::Float64, ntimes=1)
    nsys = size(Hamiltonian, 1)
    U = zeros(ComplexF64, ntimes, nsys^2, nsys^2)
    nclasstimes = trunc(Int, dt ÷ classical_dt)
    classical_dt = dt / nclasstimes
    ntotal_classical_steps = nclasstimes * ntimes
    energy = Solvents.propagate_trajectory(solvent, ps, classical_dt, ntotal_classical_steps)
    last = 0
    for t = 1:ntimes
        Utmp = exp(-1im * (Hamiltonian .+ diagm(energy[t+last,:])) * classical_dt / 2.0)
        for j = 1:nclasstimes-1
            Utmp = exp(-1im * (Hamiltonian .+ diagm(energy[t+last+j,:])) * classical_dt) * Utmp
        end
        Utmp = exp(-1im * (Hamiltonian .+ diagm(energy[t+last+nclasstimes,:])) * classical_dt / 2.0) * Utmp
        last += nclasstimes-1

        early_count = 0
        for s0p = 1:nsys
            for s0m = 1:nsys
                early_count += 1
                late_count = 0
                for s1p = 1:nsys
                    for s1m = 1:nsys
                        late_count += 1
                        U[t, late_count, early_count] += Utmp[s1p, s0p] * Utmp'[s0m, s1m]
                    end
                end
            end
        end
    end
    U
end

function calculate_average_reference_propagators(; Hamiltonian::Matrix{ComplexF64}, solvent::Solvents.Solvent, classical_dt::Float64, dt::Float64, ntimes=1)
    nsys = size(Hamiltonian, 1)
    U = zeros(ComplexF64, ntimes, nsys^2, nsys^2)
    Ucum = zeros(ComplexF64, ntimes, nsys^2, nsys^2)
    Utmp = zeros(ComplexF64, ntimes, nsys^2, nsys^2)
    for ps in solvent
        Utmp .= calculate_reference_propagators(; Hamiltonian, solvent, ps, classical_dt, dt, ntimes)
        Ucum[1, :, :] .= Utmp[1, :, :]
        for j = 2:ntimes
            Ucum[j, :, :] .= Ucum[j-1, :, :] * Utmp[j, :, :]
        end
        U .+= Ucum
    end
    U ./ solvent.num_samples
end

function calculate_bare_propagators(; Hamiltonian::Matrix{ComplexF64}, dt::Float64, ntimes=1)
    nsys = size(Hamiltonian, 1)
    U = zeros(ComplexF64, ntimes, nsys^2, nsys^2)
    for t = 1:ntimes
        Utmp = exp(-1im * Hamiltonian * dt)
        early_count = 0
        for s0p = 1:nsys
            for s0m = 1:nsys
                early_count += 1
                late_count = 0
                for s1p = 1:nsys
                    for s1m = 1:nsys
                        late_count += 1
                        U[t, late_count, early_count] += Utmp[s1p, s0p] * Utmp'[s0m, s1m]
                    end
                end
            end
        end
    end
    return U
end

function RHS(U, params, t)
    Hamil = params[1]
    for (pot, sop) in zip(params[2], params[3])
        Hamil .+= pot(t) * sop
    end
    -1im * Hamil * U
end

function calculate_bare_propagators_external_field(; Hamiltonian::Matrix{ComplexF64}, dt::Float64, ntimes=1, external_field::Vector{T}, coupling_ops::Vector{Matrix{Float64}}) where {T<:Function}
    sdim = size(Hamiltonian, 1)
    U = zeros(ComplexF64, ntimes, sdim^2, sdim^2)
    ndivs = 1000
    delt = dt / ndivs
    for t = 1:ntimes
        Utmp = Matrix{ComplexF64}(I, sdim, sdim)
        for j = 1:ndivs
            H = copy(Hamiltonian)
            for (V,s) in zip(external_field, coupling_ops)
                H .+= V(((t-1)*ndivs + (j-1)) * delt) .* s
            end
            Utmp = exp(-1im * H * delt) * Utmp
        end

        early_count = 0
        for s0p = 1:sdim
            for s0m = 1:sdim
                early_count += 1
                late_count = 0
                for s1p = 1:sdim
                    for s1m = 1:sdim
                        late_count += 1
                        U[t, late_count, early_count] += Utmp[s1p, s0p] * Utmp'[s0m, s1m]
                    end
                end
            end
        end
    end
    return U
end

end
