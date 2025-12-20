using ITensors, ITensorMPS
using LinearAlgebra
import Base

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

struct Axis{axis} end
macro Axis_str(s)
    return :(Axis{$(Expr(:quote, Symbol(s)))})
end

mutable struct MSTNPINetwork
    data::Matrix{ITensor}
end
mpoview(::Axis"time", net::MSTNPINetwork, ind::Int64) = MPO(net.data[ind, :])
mpoview(::Axis"space", net::MSTNPINetwork, ind::Int64) = MPO(net.data[:, ind])
mpsview(::Axis"time", net::MSTNPINetwork, ind::Int64) = MPS(net.data[ind, :])
mpsview(::Axis"space", net::MSTNPINetwork, ind::Int64) = MPS(net.data[:, ind])
Base.size(net::MSTNPINetwork) = size(net.data)
Base.size(net::MSTNPINetwork, j) = size(net.data, j)
Base.length(net::MSTNPINetwork) = length(net.data)

function contract_mstnpi_network(ρ0, net::MSTNPINetwork, mstnpi, tnargs::Utilities.TensorNetworkArgs)
    ρcon = apply(mpoview(Axis"time"(), net, 1), ρ0; maxdim=tnargs.maxdim, cutoff=tnargs.cutoff, alg=tnargs.algorithm)
    for j = 2:size(net, 1)-1
        ρcon = apply(mpoview(Axis"time"(), net, j), ρcon; maxdim=tnargs.maxdim, cutoff=tnargs.cutoff, alg=tnargs.algorithm)
        for k in eachindex(ρcon)
            ρcon[k] *= delta(mstnpi.sites_fb[j, k])
        end
    end
    apply(mpoview(Axis"time"(), net, size(net, 1)), ρcon; maxdim=tnargs.maxdim, cutoff=tnargs.cutoff, alg=tnargs.algorithm)
end

function init_mstnpi_network(FBProp, mstnpi, t=0)
    fbprop = deepcopy(FBProp)
    oldlinds = linkinds(fbprop)
    linds = [addtags(l, "t=$t, n=$(n)->$(n+1)") for (n, l) in enumerate(oldlinds)]
    for j in eachindex(fbprop)
        replaceinds!(fbprop[j], oldlinds, linds)
    end
    network = Matrix{ITensor}(undef, 2, mstnpi.Nsites)
    network[1,1], network[2,1] = factorize(fbprop[1], (mstnpi.sites_fb[t+1, 1], linds[1]), tags="Link,n=1,t=$t->$(t+1)", ortho="none", which_decomp="svd")
    network[1,end], network[2,end] = factorize(fbprop[end], (mstnpi.sites_fb[t+1, end], linds[end]), tags="Link,n=$(mstnpi.Nsites),t=$t->$(t+1)", ortho="none", which_decomp="svd")
    for j = 2:mstnpi.Nsites-1
        network[1,j], network[2,j] = factorize(fbprop[j], (mstnpi.sites_fb[t+1, j], linds[j-1], linds[j]), tags="Link,n=$(j),t=$t->$(t+1)", ortho="none", which_decomp="svd")
    end
    MSTNPINetwork(network)
end

function extend_mstnpi_network(netold::MSTNPINetwork, FBProp, mstnpi)
    ntimes, nsites = size(netold)
    netnew = Matrix{ITensor}(undef, ntimes+1, nsites)
    netnew[1:ntimes-1, :] .= netold.data[1:ntimes-1, :]
    nextnet = init_mstnpi_network(FBProp, mstnpi, ntimes-1)
    for j = 1:mstnpi.Nsites
        netnew[ntimes, j] = Utilities.noncontracting_prod(netold.data[ntimes, j], nextnet.data[1, j])
        netnew[ntimes+1, j] = nextnet.data[2, j]
    end
    MSTNPINetwork(netnew)
end
