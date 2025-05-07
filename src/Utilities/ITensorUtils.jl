using LinearAlgebra
using ITensors, ITensorMPS
using ITensors: scalartype

@doc raw"""
    operator_times_identity(osum::Opsum, space, iden_space, combiners)
Calculates as an MPO the the operator corresponding to ``osum \otimes 1`` where the first space is given by the indices in `space` and the second space is given by the indices in `iden_space`. Finally the MPO is put together in a single forward-backward MPO using the `combiners` .
"""
function operator_times_identity(osum::OpSum, space, iden_space, combiners)
    opmpo = MPO(osum, space) * identity_MPO(iden_space)
    for j in eachindex(opmpo)
        opmpo[j] = opmpo[j] * combiners[j]' * dag(combiners[j])
    end
    opmpo
end

"""
    calculate_Liouvillian(H::OpSum, forward_sites, backward_sites, combiners)
Returns the Liouvillian MPO corresponding to the Hamiltonian provided as an ITensor `OpSum`, to be built on the `forward_sites` and `backward_sites`, and then combined using the combiners provided.
"""
function calculate_Liouvillian(H::OpSum, forward_sites, backward_sites, combiners)
    hplus = operator_times_identity(H, forward_sites, backward_sites, combiners)
    hminus = operator_times_identity(H, backward_sites, forward_sites, combiners)
    -1im * (hplus - hminus)
end

function calculate_Liouvillian_list(H::OpSum, forward_sites, backward_sites, combiners)
    hplus = -1im * MPO(H, forward_sites) * identity_MPO(backward_sites)
    hminus = 1im * MPO(H, backward_sites) * identity_MPO(forward_sites)
    for j in eachindex(hplus)
        hplus[j] = hplus[j] * combiners[j] * combiners[j]'
        hminus[j] = hminus[j] * combiners[j] * combiners[j]'
    end
    [hplus, hminus]
end

"""
    forward_backward_combiner(sites1::Vector, sites2::Vector)
Generates combiners for generating forward-backward paths in the Liouville space from separate forward and backward path indices.
"""
function forward_backward_combiner(sites1::Vector, sites2::Vector)
    @assert length(sites1) == length(sites2) "Length of the site-sets have to be equal to generate Forward-Backward combiners."
    nsites = length(sites1)
    combiners = Vector{ITensor}(undef, nsites)
    for (n, (s1, s2)) in enumerate(zip(sites1, sites2))
        tp = tags(s1)
        tm = tags(s2)
        c = combiner(s1, s2)
        cl = combinedind(c)
        common_tags = intersect(tp, tm)
        push!(common_tags, "FBSite")
        clnew = settags(cl, join(common_tags, ", "))
        swapinds!(c, [cl], [clnew])
        combiners[n] = c
    end
    combiners
end

function forward_backward_combiner(sites1::Matrix, sites2::Matrix)
    ntimes = size(sites1, 1)
    combiners = Matrix{ITensor}(undef, size(sites1))
    combiners[1, :] = forward_backward_combiner(sites1[1, :], sites2[1, :])
    for t = 2:ntimes
        combiners[t, :] .= replacetags.(combiners[1, :], "t=0", "t=$(t-1)")
    end
    combiners
end

"""
    FBSites
Stores site information for density matrix simulations.
- `sites_forward`: Forward path sites tagged with +
- `sites_backward`: Backward path sites tagged with -
- `fb_combiners`: ITensor combiners for combining the forward and backward sites
- `sites_fb`: Forward-backward sites tagged FBSite
"""
struct FBSites
    sites_forward
    sites_backward
    fb_combiners
    sites_fb
end
"""
    fb_siteinds(sitetype, Nsites::Int; kwargs...)
Generates sites necessary for simulations in the forward-backward Liouville space. This function returns variables of the FBSites type, with forward sites being tagged with + and backward sites being tagged with -. The forward-backward sites are tagged as FBSite, and a combiner is also provided.

Arguments:
- `sitetype`: Type of the individual sites. Supports all standard ITensor site types and some custom ones defined here
- `Nsites`: Number of sites
- `kwargs...`: Other keyword arguments, most important among which are whether quantum numbers should be used
"""
function fb_siteinds(sitetype, Nsites::Int; kwargs...)
    splus = siteinds(sitetype, Nsites; addtags="+", kwargs...)
    sminus = [replacetags(dag(s), "+", "-") for s in splus]
    fbcomb = forward_backward_combiner(splus, sminus)
    FBSites(splus, sminus, fbcomb, combinedind.(fbcomb))
end

"""
    density_matrix_mps(ψ::MPS, fbsites::FBSites)
Convert a given MPS wave function defined on `fbsites.sites_forward` into a density matrix represented by an MPS over the forward-backward sites.
"""
function density_matrix_mps(ψ::MPS, fbsites::FBSites)
    ρ = MPO(ψ)
    for j in eachindex(ρ)
        ρ[j] = replaceinds(ρ[j], [fbsites.sites_forward[j]', dag(fbsites.sites_forward[j])], [fbsites.sites_forward[j], fbsites.sites_backward[j]]) * fbsites.fb_combiners[j]
    end
    convert(MPS, ρ)
end

function noncontracting_prod(A::ITensor, B::ITensor)
    Ainds = uniqueinds(A, B)
    Binds = uniqueinds(B, A)
    ABinds = commoninds(A, B)

    ainds_combiner = combiner(Ainds)
    combined_ainds = combinedind(ainds_combiner)
    binds_combiner = combiner(Binds)
    combined_binds = combinedind(binds_combiner)
    abinds_combiner = combiner(ABinds)
    combined_abinds = combinedind(abinds_combiner)

    Acomb = A * ainds_combiner * abinds_combiner
    Bcomb = B * binds_combiner * abinds_combiner
    C = ITensor(combined_ainds, combined_binds, combined_abinds)
    for aind = 1:dim(combined_ainds), bind = 1:dim(combined_binds), abind=1:dim(combined_abinds)
        C[combined_ainds=>aind, combined_binds=>bind, combined_abinds=>abind] = Acomb[combined_ainds=>aind, combined_abinds=>abind] * Bcomb[combined_binds=>bind, combined_abinds=>abind]
    end
    C * dag(ainds_combiner) * dag(binds_combiner) * dag(abinds_combiner)
end

function contract_networks(tn1::Union{AbstractVector{ITensor}, MPS, MPO}, tn2::Union{AbstractVector{ITensor}, MPS, MPO}, tnargs::Utilities.TensorNetworkArgs)
    @assert length(tn1) == length(tn2) "The length of the two networks should be equal"
    tmp = tn1[1] * tn2[1]
    left_ind = intersect(uniqueinds(tmp, tn1[2]), uniqueinds(tmp, tn2[2]))
    ans = Vector{ITensor}(undef, length(tn2))
    ans[1], tmp = factorize(tmp, left_ind; tags=tags(commonind(tn1[1], tn1[2])), cutoff=tnargs.cutoff, maxdim=tnargs.maxdim, ortho="none", which_decomp="svd")
    for j = 2:length(tn2)-1
        tmp = tmp * tn1[j] * tn2[j]
        left_ind = intersect(uniqueinds(tmp, tn1[j+1]), uniqueinds(tmp, tn2[j+1]))
        ans[j], tmp = factorize(tmp, left_ind; tags=tags(commonind(tn1[j], tn1[j+1])), cutoff=tnargs.cutoff, maxdim=tnargs.maxdim, ortho="none", which_decomp="svd")
    end
    ans[end] = tmp * tn1[end] * tn2[end]
    for j = 1:length(tn2)-1
        jinds = inds(ans[j])
        ans[j], ans[j+1] = factorize(ans[j] * ans[j+1], jinds; tags=tags(commonind(ans[j], ans[j+1])), cutoff=tnargs.cutoff, maxdim=tnargs.maxdim, ortho="none", which_decomp="svd")
    end
    ans
end

"""
    identity_MPO(sites)
Returns the identity MPO based on the given `sites`.
"""
function identity_MPO(sites)
    id = MPO(sites)
    links = linkinds(id)
    for l = 1:dim(sites[1])
        id[1][sites[1]=>l, sites[1]'=>l, links[1]=>1] = 1.0
        id[end][sites[end]=>l, sites[end]'=>l, links[end]=>1] = 1.0
    end
    for j=2:length(id)-1
        for l = 1:dim(sites[j])
            id[j][sites[j]=>l, sites[j]'=>l, links[j-1]=>1, links[j]=>1] = 1.0
        end
    end
    id
end

"""
    MPO_to_MPS(ρ::MPO, fbcombiner)
Convert a given MPO to an MPS by combining the two site indices on every tensor into a single one using the vector of combiner tensors.
"""
function MPO_to_MPS(ρ::MPO, fbcombiner)
    ρfb = deepcopy(ρ)
    for j = 1:length(ρ)
        ρfb[j] *= fbcombiner[j]
    end
    convert(MPS, ρfb)
end

"""
    MPS_to_MPO(ρ::MPS, fbcombiner)
Split a given MPS to an MPO by using the vector of combiner tensors.
"""
function MPS_to_MPO(ρ::MPS, fbcombiner)
    ρfb = deepcopy(ρ)
    for j = 1:length(ρ)
        ρfb[j] *= dag(fbcombiner[j])
    end
    convert(MPO, ρfb)
end

"""
    ITensorMPS.expect(ρ::MPO, ops; kwargs...)
Extends the ITensorMPS `expect` function to handle density matrices in the form of MPOs.
"""
function ITensorMPS.expect(ρ::MPO, ops::Tuple{<:AbstractString}; kwargs...)
    ρtmp = deepcopy(ρ)
    N = length(ρ)
    indtype = eltype(siteinds(ρ)[1])
    s = Vector{indtype}()
    sminus = Vector{indtype}()
    for ρind in siteinds(ρ)
        for j in ρind
            if hastags(j, "-")
                push!(sminus, j)
            else
                push!(s, j)
            end
        end
    end
    for j = 1:N
        swapinds!(ρtmp[j], sminus[j], dag(s[j]'))
    end

    if haskey(kwargs, :site_range)
        @warn "The `site_range` keyword arg. to `expect` is deprecated: use the keyword `sites` instead"
        sites = kwargs[:site_range]
    else
        sites = get(kwargs, :sites, 1:N)
    end

    site_range = (sites isa AbstractRange) ? sites : collect(sites)
    Ns = length(site_range)
    start_site = first(site_range)

    ElT = scalartype(ρ)
    el_types = map(o -> ishermitian(op(o, s[start_site])) ? real(ElT) : ElT, ops)

    ex = map((o, el_t) -> zeros(el_t, Ns), ops, el_types)
    for (entry, j) in enumerate(site_range)
        for (n, opname) in enumerate(ops)
            ans = if j == 1
                la = swapprime(ρtmp[j] * op(opname, s[j]), 2 => 1)
                for k = 2:N
                    la *= ρtmp[k] * delta(dag(s[k]), s[k]')
                end
                la
            else
                la = ρtmp[1] * delta(dag(s[1]), s[1]')
                for k = 2:j-1
                    la *= ρtmp[k] * delta(dag(s[k]), s[k]')
                end
                la *= swapprime(ρtmp[j] * op(opname, s[j]), 2 => 1)
                for k = j+1:N
                    la *= ρtmp[k] * delta(dag(s[k]), s[k]')
                end
                la
            end
            val = scalar(ans)
            ex[n][entry] = (el_types[n] <: Real) ? real(val) : val
        end
    end

    if sites isa Number
        return map(arr -> arr[1], ex)
    end
    return ex
end

function ITensorMPS.expect(psi::MPO, op1::AbstractString, ops::AbstractString...; kwargs...)
    return ITensorMPS.expect(psi, (op1, ops...); kwargs...)
end

function ITensorMPS.expect(psi::MPO, op1::Matrix{<:Number}, ops::Matrix{<:Number}...; kwargs...)
    return ITensorMPS.expect(psi, (op1, ops...); kwargs...)
end

function ITensorMPS.expect(psi::MPO, op::AbstractString; kwargs...)
    return first(ITensorMPS.expect(psi, (op,); kwargs...))
end

function ITensorMPS.expect(psi::MPO, op::Matrix{<:Number}; kwargs...)
    return first(ITensorMPS.expect(psi, (op,); kwargs...))
end

function build_path_amplitude_mps(fbU, sites)
    @assert length(sites)==2 "Initial path amplitude MPS should have only 2 sites corresponding to the single initial time-step."
    fbUtens = ITensor(fbU, sites[2], sites[1])
    U, V = factorize(fbUtens, (sites[1]); ortho="left", which_decomp="svd", cutoff=0.0)
    ans = MPS(2)
    ans[1] = U
    ans[2] = V
    ans
end

function extend_path_amplitude_mps(pamps, fbU, sites)
    additional_part = build_path_amplitude_mps(fbU, sites)
    ans = MPS(length(pamps) + 1)
    for (i, pa) in enumerate(pamps)
        ans[i] = deepcopy(pa)
    end
    pamps_end = deepcopy(pamps[end])
    ans[end-1] = ITensor(unioninds(pamps_end, additional_part[1]))
    old_last_link = linkinds(pamps)[end]
    site_ind = siteinds(pamps)[end]
    new_last_link = linkinds(additional_part)[1]
    for ol = 1:dim(old_last_link)
        for s = 1:dim(site_ind)
            for nl = 1:dim(new_last_link)
                ans[end-1][old_last_link=>ol, site_ind=>s, new_last_link=>nl] = pamps_end[old_last_link=>ol, site_ind=>s] * additional_part[1][site_ind=>s, new_last_link=>nl]
            end
        end
    end
    ans[end] = additional_part[2]
    ans
end

function extend_path_amplitude_mps_beyond_memory(pamps, fbU, sites)
    old_sites = siteinds(pamps)[1:3]
    trace_op = ITensor(old_sites[2])
    for iv in eachindval(old_sites[2])
        trace_op[iv] = 1
    end
    tensor13 = pamps[1] * pamps[2] * pamps[3] * trace_op
    U, V = factorize(tensor13, (old_sites[1]); ortho="left", which_decomp="svd")
    additional_part = build_path_amplitude_mps(fbU, sites)
    ans = MPS(length(pamps))
    ans[1] = U
    ans[2] = V
    for (i, pa) in enumerate(pamps[4:end])
        ans[i+2] = deepcopy(pa)
    end
    ans[end-1] = ITensor(unioninds(pamps[end], additional_part[1]))
    old_last_link = linkinds(pamps)[end]
    site_ind = siteinds(pamps)[end]
    new_last_link = linkinds(additional_part)[1]
    for ol = 1:dim(old_last_link)
        for s = 1:dim(site_ind)
            for nl = 1:dim(new_last_link)
                ans[end-1][old_last_link=>ol, site_ind=>s, new_last_link=>nl] = pamps[end][old_last_link=>ol, site_ind=>s] * additional_part[1][site_ind=>s, new_last_link=>nl]
            end
        end
    end
    ans[end] = additional_part[2]
    ans
end

function path_amplitude_to_propagator(pamps)
    ans = pamps[1]
    sinds = siteinds(pamps)
    curr_site = sinds[1]
    trace_op = ITensor(curr_site)
    for iv in eachindval(curr_site)
        trace_op[iv] = 1
    end
    for (i, ps) in enumerate(pamps[2:end-1])
        swapinds!(trace_op, [sinds[i]], [sinds[i+1]])
        ans *= ps * trace_op
    end
    noprime(ans * pamps[end])
end

function apply_contract_propagator(pamps, ifmpo)
    ans = pamps[1] * ifmpo[1]
    sinds = siteinds(pamps)
    curr_site = sinds[1]
    trace_op = ITensor(curr_site)
    for iv in eachindval(curr_site)
        trace_op[iv] = 1
    end
    for (i, ps) in enumerate(pamps[2:end-1])
        swapinds!(trace_op, [sinds[i]], [sinds[i+1]])
        ans *= ps * ifmpo[i+1] * trace_op'
    end
    noprime(ans * pamps[end] * ifmpo[end])
end

""" 
    convert_ITensor_to_matrix(tens, sinit, sterm)
Converts an ITensor with two indices to a matrix. The index `sinit` is mapped to the column and `sterm` is mapped to the row.
"""
function convert_ITensor_to_matrix(tens, sinit, sterm)
    @assert length(inds(tens)) == 2 "tens should have exactly 2 indices"
    matrix = zeros(eltype(tens), dim(sterm), dim(sinit))
    for j = 1:dim(sterm)
        for k = 1:dim(sinit)
            matrix[j, k] = tens[sinit=>k, sterm=>j]
        end
    end
    matrix
end