module TTM

using ..EtaCoefficients, ..SpectralDensities, ..Blip, ..Utilities
using ..Plot

const references = """(1) Cerrillo, J.; Cao, J. Non-Markovian Dynamical Maps: Numerical Processing of Open Quantum Trajectories. Phys. Rev. Lett. 2014, 112 (11), 110401. https://doi.org/10.1103/PhysRevLett.112.110401."""

"""
    get_propagators_QuAPI(; fbU::Array{ComplexF64,3}, Jw::Vector{T}, β, dt, ntimes, rmax, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
Calculates a timeseries of forward-backward propagators for an open quantum system using a generalized TTM fit for QuAPI. It calls the `path_integral_routine` with the bare system's forward-backward propagator and the spectral density to obtain the propagators till `rmax` time-points. Then it uses TTM to generate the other propagators.
"""
function get_propagators_QuAPI(; fbU::Array{ComplexF64,3}, Jw::Vector{T}, β, dt, ntimes, rmax, kmax::Union{Int,Nothing}=nothing, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
    @inbounds begin
        U0e_within_r, U0m_within_r, Ume_within_r, Umn_within_r = path_integral_routine(; fbU, Jw, β, dt, ntimes=rmax, kmax, extraargs, svec, verbose, reference_prop)
        T0e = zero(U0e_within_r)
        Tme = zero(Ume_within_r)
        Tmn = zero(Umn_within_r)
        for n = 1:rmax
            T0e[n, :, :] .= U0e_within_r[n, :, :]
            Tme[n, :, :] .= Ume_within_r[n, :, :]
            Tmn[n, :, :] .= Umn_within_r[n, :, :]
            for j = 1:n-1
                T0e[n, :, :] .-= Tme[j, :, :] * U0m_within_r[n-j, :, :]
                Tme[n, :, :] .-= Tme[j, :, :] * Umn_within_r[n-j, :, :]
                Tmn[n, :, :] .-= Tmn[j, :, :] * Umn_within_r[n-j, :, :]
            end
        end

        sdim2 = size(fbU, 2)
        U0e = zeros(ComplexF64, ntimes, sdim2, sdim2)
        U0m = zeros(ComplexF64, ntimes, sdim2, sdim2)
        U0e[1:rmax, :, :] = U0e_within_r
        U0m[1:rmax, :, :] = U0m_within_r
        for j = rmax+1:ntimes
            for r = 1:rmax
                U0e[j, :, :] += Tme[r, :, :] * U0m[j-r, :, :]
                U0m[j, :, :] += Tmn[r, :, :] * U0m[j-r, :, :]
            end
        end
    end
    U0e
end

function get_improved_Ts(; fbU::Array{ComplexF64,3}, Jw::Vector{T}, β, dt, ntimes, rmax, kmax::Union{Int,Nothing}=nothing, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
    ηs = [EtaCoefficients.calculate_η(jw; β, dt, kmax=ntimes, imaginary_only=reference_prop) for jw in Jw]
    _, _, _, sbar, Δs = Blip.setup_simulation(svec)
    U0e_within_r = path_integral_routine(; fbU, Jw, β, dt, ntimes=rmax, kmax, extraargs, svec, verbose, reference_prop)
    sdim2 = size(fbU, 2)
    T0e = zeros(ComplexF64, ntimes, sdim2, sdim2)
    U0e = zero(T0e)
    for n = 1:rmax
        U0e[n, :, :] .= U0e_within_r[n, :, :]
        T0e[n, :, :] .= U0e[n, :, :]
        for j = 1:n-1
            T0e[n, :, :] .-= T0e[j, :, :] * U0e[n-j, :, :]
        end
    end
    for n = rmax+1:ntimes
        U0e[n, :, :] .= sum([T0e[j, :, :] * U0e[n-j, :, :] for j = 1:n-1])
        for s1 = 1:sdim2, s2 = 1:sdim2
            val = 0
            for j = 1:length(Jw)
                val -= Δs[j, s1] * (real(ηs[j].η0e[n]) * Δs[j, s2] + 2im * imag(ηs[j].η0e[n]) * sbar[j, s2])
            end
            U0e[n, s1, s2] *= exp(val)
        end
        T0e[n, :, :] .= U0e[n, :, :]
        for j = 1:n-1
            T0e[n, :, :] .-= T0e[j, :, :] * U0e[n-j, :, :]
        end
    end
    U0e, T0e
end

function get_Ts(; fbU::Array{ComplexF64,3}, Jw::Vector{T}, β, dt, rmax, kmax::Union{Int,Nothing}=nothing, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
    U0e_within_r = path_integral_routine(; fbU, Jw, β, dt, ntimes=rmax, kmax, extraargs, svec, verbose, reference_prop)
    T0e = zero(U0e_within_r)
    for n = 1:rmax
        T0e[n, :, :] .= U0e_within_r[n, :, :]
        for j = 1:n-1
            T0e[n, :, :] .-= T0e[j, :, :] * U0e_within_r[n-j, :, :]
        end
    end
    U0e_within_r, T0e
end

function get_propagators_from_Ts(Ts, ntimes)
    ηs = [EtaCoefficients.calculate_η(jw; β, dt, kmax=ntimes, imaginary_only=reference_prop) for jw in Jw]
    _, _, _, sbar, Δs = Blip.setup_simulation(svec)
    sdim2 = size(Ts, 2)
    U0e = zeros(ComplexF64, ntimes, sdim2, sdim2)
    U0e[1, :, :] = Ts[1, :, :]
    for j = 2:ntimes
        for k = 1:min(j - 1, size(Ts, 1))
            U0e[j, :, :] += Ts[k, :, :] * U0e[j-k, :, :]
        end
    end
    U0e
end

function fit(time, vals, full_time)
    dt = time[end] - time[end-1]
    first_deriv = [-1.0 / 3, 3.0 / 2, -3.0, 11.0 / 6]
    second_deriv = [-1.0, 4.0, -5.0, 2.0]
    fhp = sum(first_deriv .* vals[end-3:end]) / dt
    fhpp = sum(second_deriv .* vals[end-3:end]) / dt^2
    α = -fhp / vals[end]
    β = -fhpp / (2 * vals[end]) + α^2 / 2
    if β < 0
        β *= α < 0 ? -0.5 : 0
    end
    index = findfirst(x -> x > time[end], full_time)
    x = full_time[index:end] .- time[end]
    argument = α * x + β * x .^ 2
    vcat(vals, vals[end] * exp.(-argument))
end

"""
    get_propagators(; fbU::Array{ComplexF64,3}, Jw::Vector{T}, β, dt, ntimes, rmax, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
Calculates a timeseries of forward-backward propagators for an open quantum system using base TTM. It calls the `path_integral_routine` with the bare system's forward-backward propagator and the spectral density to obtain the propagators till `rmax` time-points. Then it uses TTM to generate the other propagators.
"""
function get_propagators(; fbU::Array{ComplexF64,3}, Jw::Vector{T}, β, dt, ntimes, rmax, kmax::Union{Int,Nothing}=nothing, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], verbose::Bool=false, reference_prop=false, extraterm=true, fit_T=false, plot_T=false) where {T<:SpectralDensities.SpectralDensity}
    sdim2 = size(fbU, 2)
    @inbounds begin
        U0e_within_r = path_integral_routine(; fbU, Jw, β, dt, ntimes=rmax, kmax, extraargs, svec, verbose, reference_prop)
        T0e = zeros(ComplexF64, ntimes, sdim2, sdim2)
        for n = 1:rmax
            T0e[n, :, :] .= U0e_within_r[n, :, :]
            for j = 1:n-1
                T0e[n, :, :] .-= T0e[j, :, :] * U0e_within_r[n-j, :, :]
            end
        end

        U0e = zeros(ComplexF64, ntimes, sdim2, sdim2)
        U0e[1:rmax, :, :] = U0e_within_r
        t = 0:dt:(ntimes-1)*dt
        if extraterm
            ηs = [EtaCoefficients.calculate_η(jw; β, dt, kmax=ntimes, imaginary_only=reference_prop) for jw in Jw]
            _, _, _, sbar, Δs = Blip.setup_simulation(svec)
            if fit_T
                for i = 1:sdim2, j = 1:sdim2
                    T0e[:, i, j] .= fit(t[1:rmax], real.(T0e[1:rmax, i, j]), t) + 1im * fit(t[1:rmax], imag.(T0e[1:rmax, i, j]), t)
                end
                for j = rmax+1:ntimes
                    U0e[j, :, :] .= sum([T0e[r, :, :] * U0e[j-r, :, :] for r = 1:j-1])
                    for s1 = 1:sdim2, s2 = 1:sdim2
                        val = 0
                        for nb = 1:length(Jw)
                            val -= Δs[nb, s1] * (real(ηs[nb].η0e[j]) * Δs[nb, s2] + 2im * imag(ηs[nb].η0e[j]) * sbar[nb, s2])
                        end
                        U0e[j, s1, s2] *= exp(val)
                    end
                end
            else
                for j = rmax+1:ntimes
                    U0e[j, :, :] .= sum([T0e[r, :, :] * U0e[j-r, :, :] for r = 1:j-1])
                    for s1 = 1:sdim2, s2 = 1:sdim2
                        val = 0
                        for nb = 1:length(Jw)
                            val -= Δs[nb, s1] * (real(ηs[nb].η0e[j]) * Δs[nb, s2] + 2im * imag(ηs[nb].η0e[j]) * sbar[nb, s2])
                        end
                        U0e[j, s1, s2] *= exp(val)
                    end
                    T0e[j, :, :] .= U0e[j, :, :]
                    for l = 1:j-1
                        T0e[j, :, :] .-= T0e[l, :, :] * U0e[j-l, :, :]
                    end
                end
            end
        else
            if fit_T
                for i = 1:sdim2, j = 1:sdim2
                    T0e[:, i, j] .= fit(t[1:rmax], real.(T0e[1:rmax, i, j]), t) + 1im * fit(t[1:rmax], imag.(T0e[1:rmax, i, j]), t)
                end
                for j = rmax+1:ntimes
                    U0e[j, :, :] .= sum([T0e[r, :, :] * U0e[j-r, :, :] for r = 1:j-1])
                end
            else
                for j = rmax+1:ntimes
                    U0e[j, :, :] .= sum([T0e[r, :, :] * U0e[j-r, :, :] for r = 1:rmax])
                end
            end
        end
    end
    if plot_T
        for i = 1:sdim2, j = 1:sdim2
            new_figure()
            plt.plot(t[1:2rmax], real.(T0e[1:2rmax, i, j]), "k", lw=0.75, label=L"\Re")
            plt.plot(t[1:2rmax], imag.(T0e[1:2rmax, i, j]), "k--", lw=0.75, label=L"\Im")
            plt.axvline(t[rmax], c="0.75", ls="--", lw=0.75)
            plt.xlabel(L"t")
            plt.ylabel(L"T_{%$i,%$j}")
            plt.legend()
            plt.savefig("T_$(i)_$(j)_rmax$(rmax)_$(fit_T).pdf"; bbox_inches="tight")
        end
    end
    U0e
end

"""
    propagate(; fbU::Array{ComplexF64, 3}, Jw::Vector{T}, β::Real, ρ0::Matrix{ComplexF64}, dt::Real, ntimes::Int, rmax::Int, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], QuAPI::Bool=false, verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
Uses TTM to propagate the dynamics starting from `ρ0`. TTM uses propagators for different time-spans and these are calculated using the `path_integral_routine`, which returns these propagators after solving the Feynman-Vernon influence functional. If `QuAPI` is set to `false`, the default TTM method is used. Setting `QuAPI` to `true` lifts the time-translational invariance requirements of the method. Currently it is possible to use `QuAPI`, `Blip`, `TEMPO`, and `PCTNPI` to generate the propagators when `QuAPI=false`. The functions are called `build_augmented_propagators`. The additional propagators required when `QuAPI=true` can be simulated using the `build_augmented_propagators_QuAPI_TTM` of `QuAPI` and `Blip` modules.

Unlike the base methods, `TTM.propagate` cannot assume the default type of `extraargs` required for the `path_integral_routine`. Therefore, unlike `QuAPI.propagate` or `QuAPI.build_augmented_propagator`, `TTM.propagate` needs to be supplied an `extraargs` parameter compatible with the `path_integral_routine` passed in. Passing in incompatible `extraargs`, eg. `Blip.BlipArgs` with `QuAPI.build_augmented_propagator`, would result in errors.
"""
function propagate(; fbU::Array{ComplexF64,3}, Jw::Vector{T}, β::Real, ρ0::Matrix{ComplexF64}, dt::Real, ntimes::Int, rmax::Int, kmax::Union{Int,Nothing}=nothing, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], QuAPI::Bool=false, verbose::Bool=false, reference_prop=false, extraterm=true, fit_T=false, plot_T=false) where {T<:SpectralDensities.SpectralDensity}
    U0e = QuAPI ? get_propagators_QuAPI(; fbU, Jw, β, dt, ntimes, rmax, kmax, extraargs, svec, verbose, path_integral_routine, reference_prop) : get_propagators(; fbU, Jw, β, dt, ntimes, rmax, kmax, extraargs, svec, verbose, path_integral_routine, reference_prop, extraterm, fit_T, plot_T)
    Utilities.apply_propagator(; propagators=U0e, ρ0, ntimes, dt)
end

end