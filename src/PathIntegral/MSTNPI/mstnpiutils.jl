using ITensors, ITensorMPS
using LinearAlgebra

using ....Utilities

reset_time(mps::MPS; tinit=1, tnew=0) = replacetags(mps, "t=$(tinit)", "t=$(tnew)")

function get_sites(d, N, kmax, tag; kwargs...)
    tmpsite = siteinds(d, N; addtags="$(tag)", kwargs...)
    sites = Matrix{eltype(tmpsite)}(undef, kmax+1, N)
    sites[1, :] .= [replacetags(s, "Site", "t=0") for s in siteinds(d, N; addtags="$(tag)", kwargs...)]
    for t = 1:kmax
        sites[t+1, :] .= replacetags.(sites[1, :], "t=0", "t=$t")
    end
    sites
end

struct Setup{T}
    Nsites::Int64
    sites_plus::Matrix{T}
    sites_minus::Matrix{T}
    fbcombiners::Matrix{ITensor}
    sites_fb::Matrix{T}
    hamiltonian_indices::Vector{Int64}
end
function Setup(NSites, kmax, sitetype="S=1/2"; kwargs...)
    sites_plus = get_sites(sitetype, NSites, kmax, "+"; kwargs...)
    sites_minus = [replacetags(s, "+", "-") for s in sites_plus]
    combiners = Utilities.forward_backward_combiner(sites_plus, dag.(sites_minus))
    sites_fb = combinedind.(combiners)
    Setup(NSites, sites_plus, sites_minus, combiners, sites_fb, LinearIndices(sites_plus[1:2, :])[2, :])
end

function calculate_bare_propagators_single_step(; Hamiltonian, dt::AbstractFloat, ndivs=1, mstnpi::Setup, verbose::Bool=false, list::Bool=false, direct_steps::Bool=false, tnargs::Utilities.TensorNetworkArgs=Utilities.TensorNetworkArgs())
    N = mstnpi.Nsites
    ldims = [Index(1, "Link, fact") for _=1:mstnpi.Nsites-1]
    idtens = ITensor(Matrix{ComplexF64}(I, dim(mstnpi.sites_fb[1, 1]), dim(mstnpi.sites_fb[1, 1])), mstnpi.sites_fb[2, 1], dag(mstnpi.sites_fb[1, 1]))
    dyn_map0 = MPS(2*mstnpi.Nsites)
    dyn_map0[1], dyn_map0[2] = factorize(idtens, dag(mstnpi.sites_fb[1, 1]); ortho="none", which_decomp="svd")
    dyn_map0[2] *= delta(ldims[1])
    for j = 2:mstnpi.Nsites
        replaceinds!(idtens, [dag(mstnpi.sites_fb[1, j-1]), mstnpi.sites_fb[2, j-1]], [dag(mstnpi.sites_fb[1, j]), mstnpi.sites_fb[2,j]])
        dyn_map0[2j-1], dyn_map0[2j] = factorize(idtens, dag(mstnpi.sites_fb[1, j]); ortho="none", which_decomp="svd")
        dyn_map0[2j-1] *= delta(ldims[j-1])
        if j<mstnpi.Nsites
            dyn_map0[2j] *= delta(dag(ldims[j]))
        end
    end

    forward_sites = reshape(mstnpi.sites_plus[1:2, :], 2N)
    backward_sites = reshape(mstnpi.sites_minus[1:2, :], 2N)
    combiners = reshape(mstnpi.fbcombiners[1:2, :], 2N)

    verbose && @info "Calculating the Liouvillians"
    liouv, time_taken, _ = @timed if list
        Utilities.calculate_Liouvillian_list(Hamiltonian, forward_sites, backward_sites, combiners)
    else
        Utilities.calculate_Liouvillian(Hamiltonian, forward_sites, backward_sites, combiners)
    end
    if verbose
        @info "Liouvillian calculated in $(round(time_taken; digits=3)) sec."
        if !list
            ldims = linkdims(liouv)
            @info "Maximum link dimension: $(maximum(ldims)); Average link dimension: $(sum(ldims)/length(ldims))"
            if maximum(ldims) > tnargs.maxdim
                truncate!(liouv; maxdim=tnargs.maxdim, cutoff=tnargs.cutoff)
                @info "Maximum link dimension: $(maximum(ldims)); Average link dimension: $(sum(ldims)/length(ldims))"
            end
        end
        @info "Propagating TDVP"
    end

    dyn_map, time_taken, _ = @timed if direct_steps
        tdvp(liouv, dt, dyn_map0; cutoff=tnargs.cutoff, maxdim=tnargs.maxdim, nsteps=ndivs)
    else
        tdvp(liouv, dt/ndivs, dyn_map0; cutoff=tnargs.cutoff, maxdim=tnargs.maxdim, nsteps=1)
    end
    if verbose
        @info "TDVP done in $(round(time_taken; digits=3)) sec."
        ldims = linkdims(dyn_map)
        @info "Maximum link dimension of dynamical map: $(maximum(ldims))"
    end

    fbprop = MPO(N)
    for j = 1:N
        fbprop[j] = dyn_map[2j-1] * dyn_map[2j]
    end
    if verbose
        ldims = linkdims(fbprop)
        @info "Maximum link dimension of 1-step fbprop: $(maximum(ldims))"
    end
    if direct_steps
        fbprop
    else
        for j in eachindex(fbprop)
            swapinds!(fbprop[j], mstnpi.sites_fb[2, :], mstnpi.sites_fb[1, :]')
        end
        ans = fbprop
        for _ = 2:ndivs
            # ans = swapprime(ans' * fbprop, 2=>1)
            ans = apply(ans, fbprop; cutoff=tnargs.cutoff, maxdim=tnargs.maxdim, alg=tnargs.algorithm)
        end
        truncate!(ans; maxdim=tnargs.maxdim, cutoff=tnargs.cutoff)
        for j in eachindex(ans)
            swapinds!(ans[j], mstnpi.sites_fb[1, :]', mstnpi.sites_fb[2, :])
        end
        if verbose
            ldims = linkdims(ans)
            @info "Maximum link dimension of final fbprop: $(maximum(ldims))"
        end
        ans
    end
end

function calculate_bare_propagators(; Hamiltonian, dt::AbstractFloat, ndivs::Int64=1, ntimes::Int64=1, mstnpi::Setup, verbose::Bool=false, list::Bool=false, direct_steps::Bool=false, tnargs::Utilities.TensorNetworkArgs=Utilities.TensorNetworkArgs())
    fbprop = calculate_bare_propagators_single_step(; Hamiltonian, dt, ndivs, mstnpi, tnargs, list, verbose, direct_steps)
    fbprops = [deepcopy(fbprop) for _ = 1:ntimes]
    for j = 2:ntimes-1
        for k in eachindex(fbprops[j])
            replaceinds!(fbprops[j][k], mstnpi.sites_fb[1:2, k], mstnpi.sites_fb[j:j+1, k])
        end
    end
    fbprops
end