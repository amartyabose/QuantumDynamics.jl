module TTM

using ..EtaCoefficients, ..SpectralDensities, ..Utilities

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

"""
    get_propagators(; fbU::Array{ComplexF64,3}, Jw::Vector{T}, β, dt, ntimes, rmax, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
Calculates a timeseries of forward-backward propagators for an open quantum system using base TTM. It calls the `path_integral_routine` with the bare system's forward-backward propagator and the spectral density to obtain the propagators till `rmax` time-points. Then it uses TTM to generate the other propagators.
"""
function get_propagators(; fbU::Array{ComplexF64,3}, Jw::Vector{T}, β, dt, ntimes, rmax, kmax::Union{Int,Nothing}=nothing, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
    @inbounds begin
        U0e_within_r = path_integral_routine(; fbU, Jw, β, dt, ntimes=rmax, kmax, extraargs, svec, verbose, reference_prop)
        T0e = zero(U0e_within_r)
        for n = 1:rmax
            T0e[n, :, :] .= U0e_within_r[n, :, :]
            for j = 1:n-1
                T0e[n, :, :] .-= T0e[j, :, :] * U0e_within_r[n-j, :, :]
            end
        end

        sdim2 = size(fbU, 2)
        U0e = zeros(ComplexF64, ntimes, sdim2, sdim2)
        U0e[1:rmax, :, :] = U0e_within_r
        for j = rmax+1:ntimes
            for r = 1:rmax
                U0e[j, :, :] += T0e[r, :, :] * U0e[j-r, :, :]
            end
        end
    end
    U0e
end

"""
    propagate(; fbU::Array{ComplexF64, 3}, Jw::Vector{T}, β::Real, ρ0::Matrix{ComplexF64}, dt::Real, ntimes::Int, rmax::Int, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], QuAPI::Bool=false, verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
Uses TTM to propagate the dynamics starting from `ρ0`. TTM uses propagators for different time-spans and these are calculated using the `path_integral_routine`, which returns these propagators after solving the Feynman-Vernon influence functional. If `QuAPI` is set to `false`, the default TTM method is used. Setting `QuAPI` to `true` lifts the time-translational invariance requirements of the method. Currently it is possible to use `QuAPI`, `Blip`, `TEMPO`, and `PCTNPI` to generate the propagators when `QuAPI=false`. The functions are called `build_augmented_propagators`. The additional propagators required when `QuAPI=true` can be simulated using the `build_augmented_propagators_QuAPI_TTM` of `QuAPI` and `Blip` modules.

Unlike the base methods, `TTM.propagate` cannot assume the default type of `extraargs` required for the `path_integral_routine`. Therefore, unlike `QuAPI.propagate` or `QuAPI.build_augmented_propagator`, `TTM.propagate` needs to be supplied an `extraargs` parameter compatible with the `path_integral_routine` passed in. Passing in incompatible `extraargs`, eg. `Blip.BlipArgs` with `QuAPI.build_augmented_propagator`, would result in errors.
"""
function propagate(; fbU::Array{ComplexF64,3}, Jw::Vector{T}, β::Real, ρ0::Matrix{ComplexF64}, dt::Real, ntimes::Int, rmax::Int, kmax::Union{Int,Nothing}=nothing, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], QuAPI::Bool=false, verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
    U0e = QuAPI ? get_propagators_QuAPI(; fbU, Jw, β, dt, ntimes, rmax, kmax, extraargs, svec, verbose, path_integral_routine, reference_prop) : get_propagators(; fbU, Jw, β, dt, ntimes, rmax, kmax, extraargs, svec, verbose, path_integral_routine, reference_prop)
    Utilities.apply_propagator(; propagators=U0e, ρ0, ntimes, dt)
end

end