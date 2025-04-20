using ITensors, ITensorMPS
using LinearAlgebra

using ....Propagators, ....Utilities

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
    Setup(NSites, sites_plus, sites_minus, combiners, sites_fb, LinearIndices(sites_plus)[2, :])
end

function Propagators.calculate_bare_propagators(; Hamiltonian, dt::AbstractFloat, ntimes=1, mstnpi::Setup, kwargs...)
    N = mstnpi.Nsites
    idmat = Matrix{ComplexF64}(I, dim(mstnpi.sites_fb[1]), dim(mstnpi.sites_fb[1]))
    u, s, v = svd(idmat)
    links = [isodd(n) ? Index(4, "Link") : Index(1, "Link") for n = 1:2N-1]
    dyn_map0 = MPS(2N)
    dyn_map0[1] = ITensor(mstnpi.sites_fb[1], links[1])
    for s = 1:dim(mstnpi.sites_fb[1]), k = 1:dim(links[1])
        dyn_map0[1][mstnpi.sites_fb[1]=>s, links[1]=>k] = u[s, k]
    end
    dyn_map0[end] = ITensor(mstnpi.sites_fb[end], links[end])
    for s = 1:dim(mstnpi.sites_fb[end]), k = 1:dim(links[end])
        dyn_map0[end][mstnpi.sites_fb[end]=>s, links[end]=>k] = v[k, s]
    end
    for j = 2:2N-1
        dyn_map0[j] = ITensor(links[j-1], mstnpi.sites_fb[j], links[j])
        if iseven(j)
            for s = 1:dim(mstnpi.sites_fb[j]), k = 1:dim(links[j-1])
                dyn_map0[j][mstnpi.sites_fb[j]=>s, links[j-1]=>k, links[j]=>1] = v[k,s]
            end
        else
            for s = 1:dim(mstnpi.sites_fb[j]), k = 1:dim(links[j])
                dyn_map0[j][mstnpi.sites_fb[j]=>s, links[j]=>k, links[j-1]=>1] = u[s,k]
            end
        end
    end

    forward_sites = reshape(mstnpi.sites_plus[1:2, :], 2N)
    backward_sites = reshape(mstnpi.sites_minus[1:2, :], 2N)
    combiners = reshape(mstnpi.fbcombiners[1:2, :], 2N)
    dyn_map = tdvp(Utilities.calculate_Liouvillian(Hamiltonian, forward_sites, backward_sites, combiners), dt, dyn_map0; cutoff=1e-15, nsteps=ntimes);
    fbprop = MPO(N)
    for j = 1:N
        fbprop[j] = dyn_map[2j-1] * dyn_map[2j]
    end
    fbprop
end