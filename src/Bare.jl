module Bare

using DifferentialEquations

function propagate(; Hamiltonian::Matrix{ComplexF64}, ρ0::Matrix{ComplexF64}, dt::Real, ntimes::Int, cutoff=0.0, verbose::Bool=false)
    f(ρ, H, t) = -1im * (H * ρ - ρ * H')
    tspan = (0.0, ntimes * dt)
    prob = ODEProblem(f, ρ0, tspan, Hamiltonian)
    sol = solve(prob, saveat=dt)
    sdim = size(Hamiltonian, 1)
    ρs = zeros(ComplexF64, length(sol.t), sdim, sdim)
    @inbounds begin
        for j = 1:length(sol.t)
            trace = zero(ComplexF64)
            for k = 1:sdim
                trace += sol.u[j][k, k]
            end
            ρs[j, :, :] .= sol.u[j] ./ trace
        end
    end
    sol.t, ρs
end

end