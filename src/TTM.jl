module TTM

using ..EtaCoefficients, ..SpectralDensities, ..Utilities

function get_propagators(; fbU::Array{ComplexF64, 3}, Jw::Vector{T}, β, dt, ntimes, rmax, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
    @inbounds begin
        U0e_within_r, U0m_within_r, Ume_within_r, Umn_within_r = path_integral_routine(; fbU, Jw, β, dt, ntimes=rmax, extraargs, svec, verbose, reference_prop)
        T0e = similar(U0e_within_r)
        fill!(T0e, zero(ComplexF64))
        Tme = similar(Ume_within_r)
        fill!(Tme, zero(ComplexF64))
        Tmn = similar(Umn_within_r)
        fill!(Tmn, zero(ComplexF64))
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

function get_propagators_approx(; fbU::Array{ComplexF64, 3}, Jw::Vector{T}, β, dt, ntimes, rmax, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
    @inbounds begin
        U0e_within_r = path_integral_routine(; fbU, Jw, β, dt, ntimes=rmax, extraargs, svec, verbose, reference_prop, end_prop=true)
        T0e = similar(U0e_within_r)
        fill!(T0e, zero(ComplexF64))
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
                U0e[j, :, :, j] += T0e[r, :, :] * U0e[j-r, :, :]
            end
        end
    end
    U0e
end

"""
    propagate(; fbU::Array{ComplexF64, 3}, Jw::Vector{T}, β::Real, ρ0::Matrix{ComplexF64}, dt::Real, ntimes::Int, rmax::Int, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], approx::Bool=false, verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
Uses TTM to propagate the dynamics starting from `ρ0`. TTM uses propagators for different time-spans and these are calculated using the `path_integral_routine`, which returns these propagators after solving the Feynman-Vernon influence functional. Currently, one can use `QuAPI.build_augmented_propagator` and `Blip.build_augmented_propagator` with this function.
        
Unlike the base methods, `TTM.propagate` cannot assume the default type of `extraargs` required for the `path_integral_routine`. Therefore, unlike `QuAPI.propagate` or `QuAPI.build_augmented_propagator`, `TTM.propagate` needs to be supplied an `extraargs` parameter compatible with the `path_integral_routine` passed in. Passing in incompatible `extraargs`, eg. `Blip.BlipArgs` with `QuAPI.build_augmented_propagator`, would result in errors.
"""
function propagate(; fbU::Array{ComplexF64, 3}, Jw::Vector{T}, β::Real, ρ0::Matrix{ComplexF64}, dt::Real, ntimes::Int, rmax::Int, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], approx::Bool=false, verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
    U0e = approx ? get_propagators_approx(; fbU, Jw, β, dt, ntimes, rmax, extraargs, svec, verbose, path_integral_routine, reference_prop) : get_propagators(; fbU, Jw, β, dt, ntimes, rmax, extraargs, svec, verbose, path_integral_routine, reference_prop)
    Utilities.apply_propagator(; propagators=U0e, ρ0, ntimes, dt)
end

end