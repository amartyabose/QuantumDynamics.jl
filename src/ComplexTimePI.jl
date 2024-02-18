module ComplexTimePI

using LinearAlgebra
using ITensors
using FLoops
using HDF5

using ..SpectralDensities, ..Utilities, ..QuAPI, ..TEMPO

const references = """
- Topaler, M.; Makri, N. Quantum rates for a double well coupled to a dissipative bath: Accurate path integral results and comparison with approximate theories. The Journal of Chemical Physics 1994, 101 (9), 7500-7519.
- Bose, A. Quantum correlation functions through tensor network path integral. The Journal of Chemical Physics 2023, 159 (21), 214110."""

function get_time_array(t::Float64, β::Float64, N::Int64)
    Δt = t / N
    Δβ = β / 2N
    tarray = zeros(ComplexF64, 2N + 3)
    for j = 1:N
        tarray[j+1] = (j - 0.5) * (-Δt - 1im * Δβ)
    end
    tarray[N+2] = -t - 1im * β / 2
    for j = N+2:2N+1
        tarray[j+1] = (2N - j + 1.5) * (-Δt + 1im * Δβ) - 1im * β
    end
    tarray[2N+3] = -1im * β
    tarray
end

function get_B_matrix(ω::AbstractVector{<:AbstractFloat}, j::AbstractVector{<:AbstractFloat}, β::Float64, t::Float64, N::Int64)
    common_part = j ./ (ω .^ 2 .* sinh.(ω .* β ./ 2))
    tarr = get_time_array(t, β, N)
    npoints = 2N + 2
    B = zeros(ComplexF64, npoints, npoints)
    if Utilities.trapezoid(ω, common_part) ≈ 0.0
        return B
    end
    @inbounds begin
        for k = 1:npoints
            for kp = 1:k-1
                B[k, kp] = 4 / π * Utilities.trapezoid(ω, common_part .* cos.(ω .* (tarr[k+1] + tarr[k] - tarr[kp+1] - tarr[kp] + 1im * β) ./ 2) .* sin.(ω .* (tarr[k+1] - tarr[k]) ./ 2) .* sin.(ω .* (tarr[kp+1] - tarr[kp]) ./ 2))
                B[kp, k] = B[k, kp]
            end
            B[k, k] = 2 / π * Utilities.trapezoid(ω, common_part .* sin.(ω .* (tarr[k+1] - tarr[k] + 1im * β) ./ 2) .* sin.(ω .* (tarr[k+1] - tarr[k]) ./ 2))
        end
    end
    B
end

function get_B_matrix(J::SpectralDensities.SpectralDensity, β::Float64, t::Float64, N::Int64)
    ω, j = SpectralDensities.tabulate(J, false)
    get_B_matrix(ω, j, β, t, N)
end

function get_complex_time_propagator(Hamiltonian::Matrix{ComplexF64}, β::Float64, t::Float64, N::Int64)
    Δtc = (t - 1im * β / 2) / N
    exp(-1im * Hamiltonian * Δtc), exp(1im * Hamiltonian * conj(Δtc))
end

"""
    unnormalized_correlation_function_quapi(; Hamiltonian::AbstractMatrix{ComplexF64}, β::Real, t::Real, N::Int, Jw::AbstractVector{<:SpectralDensities.SpectralDensity}, svec::AbstractMatrix{<:Real}, A, B, extraargs::QuAPI.QuAPIArgs=QuAPI.QuAPIArgs())
Calculates the ``<A(0) B(t)>`` correlation function for a system interacting with an environment at a time-point `t` using QuAPI.

Relevant references:
- Topaler, M.; Makri, N. Quantum rates for a double well coupled to a dissipative bath: Accurate path integral results and comparison with approximate theories. The Journal of Chemical Physics 1994, 101 (9), 7500-7519.

Arguments:
- `Hamiltonian`: system Hamiltonian
- `Jw`: array of spectral densities
- `svec`: diagonal elements of system operators through which the corresponding baths interact. QuAPI currently only works for baths with diagonal coupling to the system.
- `β`: inverse temperature
- `t`: time at which the correlation function is evaluated
- `N`: number of path integral discretizations
- `A`: system operator evaluated at time zero
- `B`: system operator evaluated at time `t`
- `extraargs`: extra arguments for the QuAPI algorithm. Contains a `cutoff` threshold for filtration of paths.
"""
function unnormalized_correlation_function_quapi(; Hamiltonian::Matrix{ComplexF64}, β::Float64, t::Float64, N::Int64, Jw::Vector{<:SpectralDensities.SpectralDensity}, svec::Matrix{Float64}, A, B, extraargs::QuAPI.QuAPIArgs=QuAPI.QuAPIArgs())
    @assert length(Jw) == size(svec, 1)
    nbaths = length(Jw)
    U, Udag = get_complex_time_propagator(Hamiltonian, β, t, N)
    Bmat = [get_B_matrix(J.ω, J.jw, β, t, N) for J in Jw]
    npoints = size(Bmat[1], 1)
    nfor = npoints ÷ 2
    sdim = size(Hamiltonian, 1)
    num_paths = sdim^npoints
    amp = 0.0 + 0.0im
    for path_num = 1:num_paths
        states = Utilities.unhash_path(path_num, npoints - 1, sdim)
        # e^{-i H tc} A e^{i H tc} B
        # A e^{i H tc} B e^{-i H tc}
        val = A[states[nfor], states[nfor+1]] * B[states[end], states[1]]
        if abs(val) ≤ extraargs.cutoff
            continue
        end
        for j = 1:nfor-1
            val *= U[states[j], states[j+1]]
        end
        for j = nfor+1:npoints-1
            val *= Udag[states[j], states[j+1]]
        end
        if abs(val) ≤ extraargs.cutoff
            continue
        end
        infl = 0.0 + 0.0im
        for nb = 1:nbaths
            nnonzeros = 0
            for s in states
                if svec[nb, s] != 0
                    nnonzeros += 1
                end
            end
            ks = zeros(Int64, nnonzeros)
            svecs = zeros(nnonzeros)
            l = 1
            for (k, s) in enumerate(states)
                if svec[nb, s] != 0
                    svecs[l] = svec[nb, s]
                    ks[l] = k
                    l += 1
                end
            end
            for (i, k) in enumerate(ks)
                for kp = 1:i
                    infl -= Bmat[nb][k, ks[kp]] * svecs[i] * svecs[kp]
                end
            end
        end
        amp += val * exp(infl)
    end
    amp
end

function correlation_function_quapi_parallel(; Hamiltonian::Matrix{ComplexF64}, β::Float64, t::Float64, N::Int64, Jw::Vector{<:SpectralDensities.SpectralDensity}, svec::Matrix{Float64}, A, B, extraargs::QuAPI.QuAPIArgs=QuAPI.QuAPIArgs())
    @assert length(Jw) == size(svec, 1)
    nbaths = length(Jw)
    U, Udag = get_complex_time_propagator(Hamiltonian, β, t, N)
    Bmat = [get_B_matrix(J.ω, J.jw, β, t, N) for J in Jw]
    npoints = size(Bmat[1], 1)
    nfor = npoints ÷ 2
    sdim = size(Hamiltonian, 1)
    num_paths = sdim^npoints
    @floop for path_num = 1:num_paths
        states = Utilities.unhash_path(path_num, npoints - 1, sdim)
        # e^{-i H tc} A e^{i H tc} B
        # A e^{i H tc} B e^{-i H tc}
        val = A[states[nfor], states[nfor+1]] * B[states[end], states[1]]
        if abs(val) ≤ extraargs.cutoff
            continue
        end
        for j = 1:nfor-1
            val *= U[states[j], states[j+1]]
        end
        for j = nfor+1:npoints-1
            val *= Udag[states[j], states[j+1]]
        end
        if abs(val) ≤ extraargs.cutoff
            continue
        end
        infl = 0.0 + 0.0im
        for nb = 1:nbaths
            nnonzeros = 0
            for s in states
                if svec[nb, s] != 0
                    nnonzeros += 1
                end
            end
            ks = zeros(Int64, nnonzeros)
            svecs = zeros(nnonzeros)
            l = 1
            for (k, s) in enumerate(states)
                if svec[nb, s] != 0
                    svecs[l] = svec[nb, s]
                    ks[l] = k
                    l += 1
                end
            end
            for (i, k) in enumerate(ks)
                for kp = 1:i
                    infl -= Bmat[nb][k, ks[kp]] * svecs[i] * svecs[kp]
                end
            end
        end
        @reduce amp = 0.0 + val * exp(infl)
    end
    amp
end

function get_Bmat_MPO_left(svec, sites, k, Bmat, npoints)
    nsites = length(sites)
    sitedim = dim(sites[1])
    nbaths = size(svec, 1)
    nstates = size(svec, 2)
    links = [Index(nstates, "Link") for j in 1:nsites-1]
    Bmat_MPO = MPO(nsites)
    loc = nsites ÷ 2

    @inbounds begin
        tensor = ITensor(sites[1], sites[1]', links[1])
        for s = 1:sitedim, β = 1:nstates
            tensor[sites[1]=>s, sites[1]'=>s, links[1]=>β] = exp(sum([-Bmat[nb][1, k] * svec[nb, s] * svec[nb, β] for nb = 1:nbaths]))
        end
        Bmat_MPO[1] = tensor

        for (j, site) in enumerate(sites[2:loc-1])
            tensor = ITensor(links[j], site, site', links[j+1])
            for s = 1:sitedim, β = 1:nstates
                tensor[links[j]=>β, site=>s, site'=>s, links[j+1]=>β] = exp(sum([-Bmat[nb][j+1, k] * svec[nb, s] * svec[nb, β] for nb = 1:nbaths]))
            end
            Bmat_MPO[j+1] = tensor
        end

        tensor = ITensor(links[loc-1], sites[loc], links[loc])
        for s = 1:sitedim
            tensor[links[loc-1]=>s, sites[loc]=>s, links[loc]=>s] = exp(sum([-Bmat[nb][k, k] * svec[nb, s] * svec[nb, s] for nb = 1:nbaths]))
        end
        Bmat_MPO[loc] = tensor

        for (j, site) in enumerate(sites[loc+1:end-1])
            tensor = ITensor(links[j+loc-1], site, site', links[j+loc])
            for s = 1:sitedim, β = 1:nstates
                tensor[links[j+loc-1]=>β, site=>s, site'=>s, links[j+loc]=>β] = exp(sum([-Bmat[nb][k, k+j+npoints-nsites] * svec[nb, s] * svec[nb, β] for nb = 1:nbaths]))
            end
            Bmat_MPO[j+loc] = tensor
        end

        tensor = ITensor(links[end], sites[end], sites[end]')
        for s = 1:sitedim, β = 1:nstates
            tensor[links[end]=>β, sites[end]=>s, sites[end]'=>s] = exp(sum([-Bmat[nb][k, npoints] * svec[nb, s] * svec[nb, β] for nb = 1:nbaths]))
        end
        Bmat_MPO[end] = tensor
    end

    Bmat_MPO
end

function get_Bmat_MPO_right(svec, sites, k, Bmat, npoints)
    nsites = length(sites)
    sitedim = dim(sites[1])
    nbaths = size(svec, 1)
    nstates = size(svec, 2)
    links = [Index(nstates, "Link") for j in 1:nsites-1]
    Bmat_MPO = MPO(nsites)
    loc = nsites ÷ 2 + 1

    @inbounds begin
        tensor = ITensor(sites[1], sites[1]', links[1])
        for s = 1:sitedim, β = 1:nstates
            tensor[sites[1]=>s, sites[1]'=>s, links[1]=>β] = exp(sum([-Bmat[nb][1, k] * svec[nb, s] * svec[nb, β] for nb = 1:nbaths]))
        end
        Bmat_MPO[1] = tensor

        for (j, site) in enumerate(sites[2:loc-1])
            tensor = ITensor(links[j], site, site', links[j+1])
            for s = 1:sitedim, β = 1:nstates
                tensor[links[j]=>β, site=>s, site'=>s, links[j+1]=>β] = exp(sum([-Bmat[nb][j+1, k] * svec[nb, s] * svec[nb, β] for nb = 1:nbaths]))
            end
            Bmat_MPO[j+1] = tensor
        end

        tensor = ITensor(links[loc-1], sites[loc], links[loc])
        for s = 1:sitedim
            tensor[links[loc-1]=>s, sites[loc]=>s, links[loc]=>s] = exp(sum([-Bmat[nb][k, k] * svec[nb, s] * svec[nb, s] for nb = 1:nbaths]))
        end
        Bmat_MPO[loc] = tensor

        for (j, site) in enumerate(sites[loc+1:end-1])
            tensor = ITensor(links[j+loc-1], site, site', links[j+loc])
            for s = 1:sitedim, β = 1:nstates
                tensor[links[j+loc-1]=>β, site=>s, site'=>s, links[j+loc]=>β] = exp(sum([-Bmat[nb][k, k+j] * svec[nb, s] * svec[nb, β] for nb = 1:nbaths]))
            end
            Bmat_MPO[j+loc] = tensor
        end

        tensor = ITensor(links[end], sites[end], sites[end]')
        for s = 1:sitedim, β = 1:nstates
            tensor[links[end]=>β, sites[end]=>s, sites[end]'=>s] = exp(sum([-Bmat[nb][k, npoints] * svec[nb, s] * svec[nb, β] for nb = 1:nbaths]))
        end
        Bmat_MPO[end] = tensor
    end

    Bmat_MPO
end

"""
    A_of_t(; Hamiltonian::AbstractMatrix{ComplexF64}, β::Real, t::Real, N::Int, Jw::AbstractVector{<:SpectralDensities.SpectralDensity}, svec::AbstractMatrix{<:Real}, A, B, extraargs::Utilities.TensorNetworkArgs=Utilities.TensorNetworkArgs())
Calculates ``tr_env(U(t) exp(-β H/2) A exp(-β H/2) U^{-1}(t))`` for a system interacting with an environment at a time-point `t` using the tensor network path integral method. This can be used for thermodynamics or for calculating correlation functions.

Relevant references:
$(references)

Arguments:
- `Hamiltonian`: system Hamiltonian
- `Jw`: array of spectral densities
- `svec`: diagonal elements of system operators through which the corresponding baths interact. QuAPI currently only works for baths with diagonal coupling to the system.
- `β`: inverse temperature
- `t`: time at which the function is evaluated
- `N`: number of path integral discretizations
- `A`: system operator to be evaluated
- `extraargs`: extra arguments for the tensor network algorithm. Contains the `cutoff` threshold for SVD filtration, the maximum bond dimension, `maxdim`, and the `algorithm` of applying an MPO to an MPS.
"""
function A_of_t(; Hamiltonian::AbstractMatrix{ComplexF64}, β::Real, t::Real, N::Int, Jw::AbstractVector{<:SpectralDensities.SpectralDensity}, svec::AbstractMatrix{<:Real}, A, extraargs::Utilities.TensorNetworkArgs=Utilities.TensorNetworkArgs(), verbose::Bool=false)
    @assert length(Jw) == size(svec, 1)
    nbaths = length(Jw)
    U, Udag = get_complex_time_propagator(Hamiltonian, β, t, N)
    Bmat = [get_B_matrix(J, β, t, N) for J in Jw]
    npoints = size(Bmat[1], 1)
    nfor = npoints ÷ 2
    sdim = size(Hamiltonian, 1)

    sites = siteinds(sdim, npoints)
    Ufor, S, V = svd(U)
    Rfor = diagm(S) * V'

    UB, S, V = svd(A)
    RB = diagm(S) * V'

    Uback, S, V = svd(Udag)
    Rback = diagm(S) * V'

    # e^{-i H tc} A e^{i H tc} B
    pathmps = MPS(npoints)
    linkleft = Index(size(Rfor, 1), "Link")
    pathmps[1] = ITensor(Ufor, sites[1], linkleft)
    for s = 1:dim(sites[1])
        selfinfl = exp(sum([-Bmat[nb][1, 1] * svec[nb, s] * svec[nb, s] for nb = 1:nbaths]))
        for l = 1:dim(linkleft)
            pathmps[1][sites[1]=>s, linkleft=>l] *= selfinfl
        end
    end
    for j = 2:nfor-1
        linkright = Index(size(Ufor, 1), "Link")
        pathmps[j] = ITensor(linkleft, sites[j], linkright)
        for ll = 1:dim(linkleft), sl = 1:dim(sites[j]), lr = 1:dim(linkright)
            pathmps[j][linkleft=>ll, sites[j]=>sl, linkright=>lr] = Rfor[ll, sl] * Ufor[sl, lr]
        end
        linkleft = linkright
    end
    linkright = Index(size(RB, 1), "Link")
    pathmps[nfor] = ITensor(linkleft, sites[nfor], linkright)
    for ll = 1:dim(linkleft), sl = 1:dim(sites[nfor]), lr = 1:dim(linkright)
        pathmps[nfor][linkleft=>ll, sites[nfor]=>sl, linkright=>lr] = Rfor[ll, sl] * UB[sl, lr]
    end
    linkleft = linkright
    linkright = Index(size(Uback, 1), "Link")
    pathmps[nfor+1] = ITensor(linkleft, sites[nfor+1], linkright)
    for ll = 1:dim(linkleft), sl = 1:dim(sites[nfor+1]), lr = 1:dim(linkright)
        pathmps[nfor+1][linkleft=>ll, sites[nfor+1]=>sl, linkright=>lr] = RB[ll, sl] * Uback[sl, lr]
    end
    linkleft = linkright
    for j = nfor+2:npoints-1
        linkright = Index(size(Uback, 1), "Link")
        pathmps[j] = ITensor(linkleft, sites[j], linkright)
        for ll = 1:dim(linkleft), sl = 1:dim(sites[j]), lr = 1:dim(linkright)
            pathmps[j][linkleft=>ll, sites[j]=>sl, linkright=>lr] = Rback[ll, sl] * Uback[sl, lr]
        end
        linkleft = linkright
    end
    pathmps[npoints] = ITensor(Rback, linkleft, sites[npoints])
    for s = 1:dim(sites[end])
        selfinfl = exp(sum([-Bmat[nb][end, end] * svec[nb, s] * svec[nb, s] for nb = 1:nbaths]))
        for l = 1:dim(linkleft)
            pathmps[end][sites[end]=>s, linkleft=>l] *= selfinfl
        end
    end

    maxlinkdims = []
    for k = npoints÷2:-1:2
        Bmat_MPO = get_Bmat_MPO_left(svec, siteinds(pathmps), k, Bmat, npoints)
        new_pathmps = apply(Bmat_MPO, pathmps; cutoff=extraargs.cutoff, maxdim=extraargs.maxdim, alg=extraargs.algorithm)
        pathmps = MPS(length(new_pathmps) - 1)
        for s = 1:k-1
            pathmps[s] = new_pathmps[s]
        end
        pathmps[k] = new_pathmps[k] * new_pathmps[k+1]
        for s = k+1:length(pathmps)
            pathmps[s] = new_pathmps[s+1]
        end
        append!(maxlinkdims, maxlinkdim(pathmps))

        Bmat_MPO = get_Bmat_MPO_right(svec, siteinds(pathmps), npoints - k + 1, Bmat, npoints)
        new_pathmps = apply(Bmat_MPO, pathmps; cutoff=extraargs.cutoff, maxdim=extraargs.maxdim, alg=extraargs.algorithm)
        pathmps = MPS(length(new_pathmps) - 1)
        for s = 1:k-1
            pathmps[s] = new_pathmps[s]
        end
        pathmps[k] = new_pathmps[k] * new_pathmps[k+1]
        for s = k+1:length(pathmps)
            pathmps[s] = new_pathmps[s+1]
        end
        append!(maxlinkdims, maxlinkdim(pathmps))
    end

    avg_bond = sum(maxlinkdims) / length(maxlinkdims)
    if verbose
        @info "Average bond dimension = $(verbose)"
    end

    tempmat = Utilities.convert_ITensor_to_matrix(TEMPO.path_amplitude_to_propagator(pathmps), sites[end], sites[1])
    for sl = 1:sdim, sr = 1:sdim
        tempmat[sl, sr] *= exp(sum([-Bmat[nb][1, end] * svec[nb, sl] * svec[nb, sr] for nb = 1:nbaths]))
    end
    avg_bond, tempmat
end

"""
    unnormalized_correlation_function_tnpi(; Hamiltonian::AbstractMatrix{ComplexF64}, β::Real, t::Real, N::Int, Jw::AbstractVector{<:SpectralDensities.SpectralDensity}, svec::AbstractMatrix{<:Real}, A, B, extraargs::Utilities.TensorNetworkArgs=Utilities.TensorNetworkArgs())
Calculates the ``<A(0) B(t)>`` correlation function for a system interacting with an environment at a time-point `t` using the tensor network path integral method.

Relevant references:
$(references)

Arguments:
- `Hamiltonian`: system Hamiltonian
- `Jw`: array of spectral densities
- `svec`: diagonal elements of system operators through which the corresponding baths interact. QuAPI currently only works for baths with diagonal coupling to the system.
- `β`: inverse temperature
- `t`: time at which the correlation function is evaluated
- `N`: number of path integral discretizations
- `A`: system operator evaluated at time zero
- `B`: array of system operators evaluated at time `t`
- `extraargs`: extra arguments for the tensor network algorithm. Contains the `cutoff` threshold for SVD filtration, the maximum bond dimension, `maxdim`, and the `algorithm` of applying an MPO to an MPS.
"""
function single_point_correlation_function_tnpi(; Hamiltonian::AbstractMatrix{ComplexF64}, β::Real, t::Real, N::Int, Jw::AbstractVector{<:SpectralDensities.SpectralDensity}, svec::AbstractMatrix{<:Real}, A, B, extraargs::Utilities.TensorNetworkArgs=Utilities.TensorNetworkArgs(), verbose::Bool=false)
    avg_bond_dim, At = A_of_t(; Hamiltonian, β, t, N, Jw, svec, A, extraargs, verbose)
    length(B) == 1 ? (avg_bond_dim, tr(B[1] * At)) : (avg_bond_dim, [tr(b * At) for b in B])
end

"""
    correlation_function_tnpi(; Hamiltonian::AbstractMatrix{ComplexF64}, β::Real, tfinal::Real, dt::Real, N::Int, Jw::AbstractVector{<:SpectralDensities.SpectralDensity}, svec::AbstractMatrix{<:Real}, A, B, Z::Real, extraargs::Utilities.TensorNetworkArgs=Utilities.TensorNetworkArgs(), verbose::Bool=false, output::Union{Nothing,HDF5.Group}=nothing)
Calculates the ``\frac{<A(0) B(t)>}{Z}`` correlation function for a system interacting with an environment upto a maximum time of `tfinal` with a time-step of `dt` using the tensor network path integral method.

Relevant references:
$(references)

Arguments:
- `Hamiltonian`: system Hamiltonian
- `Jw`: array of spectral densities
- `svec`: diagonal elements of system operators through which the corresponding baths interact. QuAPI currently only works for baths with diagonal coupling to the system.
- `β`: inverse temperature
- `tfinal`: maximum time till which the correlation function is evaluated
- `dt`: time-step for evaluating correlation function
- `N`: number of path integral discretizations
- `A`: system operator evaluated at time zero
- `B`: array of system operators evaluated at time `t`
- `Z`: partition function for normalization
- `extraargs`: extra arguments for the tensor network algorithm. Contains the `cutoff` threshold for SVD filtration, the maximum bond dimension, `maxdim`, and the `algorithm` of applying an MPO to an MPS.
- `verbose`: verbosity
- `output`: output HDF5 file for storage of results
"""
function correlation_function_tnpi(; Hamiltonian::AbstractMatrix{ComplexF64}, β::Real, tfinal::Real, dt::Real, N::Int, Jw::AbstractVector{<:SpectralDensities.SpectralDensity}, svec::AbstractMatrix{<:Real}, A, B, Z::Real, extraargs::Utilities.TensorNetworkArgs=Utilities.TensorNetworkArgs(), verbose::Bool=false, output::Union{Nothing,HDF5.Group}=nothing)
    time = 0:dt:tfinal
    corr = zeros(ComplexF64, length(time), length(B))
    if !isnothing(output)
        Utilities.check_or_insert_value(output, "time", time)
        Utilities.check_or_insert_value(output, "corr", corr)
        Utilities.check_or_insert_value(output, "bond_dims", zeros(Float64, length(time)))
    end
    for (i, t) in enumerate(time)
        _, time_taken, memory_allocated, gc_time, _ = @timed begin
            avg_bond_dim, At = A_of_t(; Hamiltonian, β, t, N, Jw, svec, A, extraargs, verbose)
        end
        At /= Z
        corr[i, :] .= [tr(b * At) for b in B]
        if verbose
            @info "Step = $(i); avg bond dimension = $(avg_bond_dim); time = $(round(time_taken; digits=3)) sec; memory allocated = $(round(memory_allocated / 1e6; digits=3)) GB; gc time = $(round(gc_time; digits=3)) sec"
        end
        if !isnothing(output)
            output["corr"][i, :] .= corr[i, :]
            output["bond_dims"][i] = avg_bond_dim
        end
    end
    corr
end

end
