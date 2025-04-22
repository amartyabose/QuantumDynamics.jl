using ITensors, ITensorMPS
using LinearAlgebra

using ....Utilities

reset_time(mps::MPS; tinit=1, tnew=0) = replacetags(mps, "t=$(tinit)", "t=$(tnew)")

function get_sites(d, N, kmax, tag)
    sites = Matrix{Index{Int64}}(undef, kmax+1, N)
    sites[1, :] .= [replacetags(s, "Site", "t=0") for s in siteinds(d, N; addtags="$(tag)")]
    for t = 1:kmax
        sites[t+1, :] .= replacetags.(sites[1, :], "t=0", "t=$t")
    end
    sites
end

struct Setup
    Nsites::Int64
    sites_plus::Matrix{Index{Int64}}
    sites_minus::Matrix{Index{Int64}}
    fbcombiners::Matrix{ITensor}
    sites_fb::Matrix{Index{Int64}}
    hamiltonian_indices::Vector{Int64}
end
function Setup(NSites, kmax, sitetype="S=1/2")
    sites_plus = get_sites(sitetype, NSites, kmax, "+")
    sites_minus = [replacetags(s, "+", "-") for s in sites_plus]
    combiners = Utilities.forward_backward_combiner(sites_plus, sites_minus)
    sites_fb = combinedind.(combiners)
    Setup(NSites, sites_plus, sites_minus, combiners, sites_fb, LinearIndices(sites_plus[1:2, :])[2, :])
end

function calculate_bare_propagators(; Hamiltonian, dt::AbstractFloat, ntimes=1, mstnpi::Setup, verbose::Bool=false, list::Bool=false, direct_steps::Bool=false, tnargs::Utilities.TensorNetworkArgs=Utilities.TensorNetworkArgs())
    N = mstnpi.Nsites
    idmat = Matrix{ComplexF64}(I, dim(mstnpi.sites_fb[1]), dim(mstnpi.sites_fb[1]))
    u, s, v = svd(idmat)
    links = [isodd(n) ? Index(4, "Link") : Index(1, "Link") for n = 1:2N-1]
    dyn_map0 = MPS(2N)
    dyn_map0[1] = ITensor(mstnpi.sites_fb[1,1], links[1])
    for s = 1:dim(mstnpi.sites_fb[1,1]), k = 1:dim(links[1])
        dyn_map0[1][mstnpi.sites_fb[1,1]=>s, links[1]=>k] = u[s, k]
    end
    dyn_map0[end] = ITensor(mstnpi.sites_fb[2,end], links[end])
    for s = 1:dim(mstnpi.sites_fb[2,end]), k = 1:dim(links[end])
        dyn_map0[end][mstnpi.sites_fb[2,end]=>s, links[end]=>k] = v[k, s]
    end
    for j = 2:2N-1
        n = trunc(Int64, ceil(j/2))
        t = mod1(j, 2)
        dyn_map0[j] = ITensor(links[j-1], mstnpi.sites_fb[t, n], links[j])
        if iseven(j)
            for s = 1:dim(mstnpi.sites_fb[t, n]), k = 1:dim(links[j-1])
                dyn_map0[j][mstnpi.sites_fb[t, n]=>s, links[j-1]=>k, links[j]=>1] = v[k,s]
            end
        else
            for s = 1:dim(mstnpi.sites_fb[t, n]), k = 1:dim(links[j])
                dyn_map0[j][mstnpi.sites_fb[t, n]=>s, links[j]=>k, links[j-1]=>1] = u[s,k]
            end
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
        tdvp(liouv, dt, dyn_map0; cutoff=tnargs.cutoff, maxdim=tnargs.maxdim, nsteps=ntimes)
    else
        tdvp(liouv, dt/ntimes, dyn_map0; cutoff=tnargs.cutoff, maxdim=tnargs.maxdim, nsteps=1)
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
        for _ = 2:ntimes
            # ans = swapprime(ans' * fbprop, 2=>1)
            ans = apply(ans, fbprop; cutoff=tnargs.cutoff, maxdim=tnargs.maxdim, alg=tnargs.algorithm)
        end
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