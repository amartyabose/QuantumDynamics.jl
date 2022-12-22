module TTM

using ..EtaCoefficients, ..SpectralDensities

@inbounds function get_propagators(; fbU::Matrix{ComplexF64}, Jw::Vector{T}, β, dt, ntimes, rmax, build_propagator, cutoff=-1, svec=[1.0 -1.0], verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
    U0e_within_r, U0m_within_r, Ume_within_r, Umn_within_r = build_propagator(; fbU, Jw, β, dt, ntimes=rmax, cutoff, svec, verbose, reference_prop)
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

    sdim2 = size(fbU, 1)
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
    return U0e
end

@inbounds function get_propagators_approx(; fbU::Matrix{ComplexF64}, Jw::Vector{T}, β, dt, ntimes, rmax, build_propagator, cutoff=-1, svec=[1.0 -1.0], verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
    U0e_within_r, _, _, _ = build_propagator(; fbU, Jw, β, dt, ntimes=rmax, cutoff, svec, verbose, reference_prop)
    T0e = similar(U0e_within_r)
    fill!(T0e, zero(ComplexF64))
    for n = 1:rmax
        T0e[n, :, :] .= U0e_within_r[n, :, :]
        for j = 1:n-1
            T0e[n, :, :] .-= T0e[j, :, :] * U0e_within_r[n-j, :, :]
        end
    end

    sdim2 = size(fbU, 1)
    U0e = zeros(ComplexF64, ntimes, sdim2, sdim2)
    U0e[1:rmax, :, :] = U0e_within_r
    for j = rmax+1:ntimes
        for r = 1:rmax
            U0e[j, :, :, j] += T0e[r, :, :] * U0e[j-r, :, :]
        end
    end
    return U0e
end

@inbounds function propagate_using_propagators(; propagators, ρ0, ntimes)
    sdim = size(ρ0, 1)
    ρs = zeros(ComplexF64, ntimes + 1, sdim, sdim)
    ρs[1, :, :] = ρ0
    ρvec = collect(Iterators.flatten(ρ0))
    for j = 1:ntimes
        ρs[j+1, :, :] = reshape(propagators[j, :, :] * ρvec, (sdim, sdim))
    end
    ρs
end

@inbounds function propagate(; fbU::Matrix{ComplexF64}, Jw::Vector{T}, β::Real, ρ0::Matrix{ComplexF64}, dt::Real, ntimes::Int, rmax::Int, build_propagator, cutoff=-1, svec=[1.0 -1.0], approx::Bool=false, verbose::Bool=false, reference_prop=false) where {T<:SpectralDensities.SpectralDensity}
    U0e = approx ? get_propagators_approx(; fbU, Jw, β, dt, ntimes, rmax, cutoff, svec, verbose, build_propagator, reference_prop) : get_propagators(; fbU, Jw, β, dt, ntimes, rmax, cutoff, svec, verbose, build_propagator, reference_prop)
    0:dt:ntimes*dt, propagate_using_propagators(; propagators=U0e, ρ0, ntimes)
end

end