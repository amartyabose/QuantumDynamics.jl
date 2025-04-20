using LinearAlgebra
using ITensors, ITensorMPS
using ITensors: scalartype

"""
    calculate_Liouvillian(H::OpSum, sites)
Returns the forward-backward space combiner and the Liouvillian MPO corresponding to the Hamiltonian provided as an ITensor `OpSum`, to be built on the `sites`.
"""
function calculate_Liouvillian(H::OpSum, forward_sites, backward_sites, combiners)
    hplus = MPO(H, forward_sites) * identity_MPO(backward_sites)
    hminus = MPO(H, backward_sites) * identity_MPO(forward_sites)
    liouv = -1im * (hplus - hminus)
    for j in eachindex(liouv)
        liouv[j] = liouv[j] * combiners[j] * combiners[j]'
    end
    liouv
end

function forward_backward_combiner(sites1, sites2)
    nsites = size(sites1, 2)
    ntimes = size(sites1, 1)
    combiners = Matrix{ITensor}(undef, size(sites1))

    for n = 1:nsites
        tp = tags(sites1[1, n])
        tm = tags(sites2[1, n])
        c = combiner(sites1[1, n], sites2[1, n])
        cl = combinedind(c)
        common_tags = intersect(tp, tm)
        push!(common_tags, "FBSite")
        clnew = settags(cl, join(common_tags, ", "))
        swapinds!(c, [cl], [clnew])
        combiners[1, n] = c
    end
    for t = 2:ntimes
        combiners[t, :] .= replacetags.(combiners[1, :], "t=0", "t=$(t-1)")
    end
    combiners
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
    s = Vector{Index{Int64}}()
    sminus = Vector{Index{Int64}}()
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
        swapinds!(ρtmp[j], sminus[j], s[j]')
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
            ans = j == 1 ? swapprime(ρtmp[j] * op(opname, s[j])', 2 => 1) * delta(s[j], s[j]') : ρtmp[1] * delta(s[1], s[1]')
            for k = 2:j-1
                ans *= ρtmp[k] * delta(s[k], s[k]')
            end
            if j != 1
                ans *= swapprime(ρtmp[j] * op(opname, s[j])', 2 => 1) * delta(s[j], s[j]')
            end
            for k = j+1:N
                ans *= ρtmp[k] * delta(s[k], s[k]')
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
    fbUtens = ITensor(fbU, sites)
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
    matrix = zeros(eltype(tens), dim(sterm), dim(sinit))
    for j = 1:dim(sterm)
        for k = 1:dim(sinit)
            matrix[j, k] = tens[sinit=>k, sterm=>j]
        end
    end
    matrix
end

function ITensors.:+(::ITensors.Algorithm"approx", ψ::MPST...; cutoff::Float64, maxdim::Int64) where {MPST<:ITensorMPS.AbstractMPS}
    @error "This approximate function needs to be implemented"
    n = length(first(ψ))
    @assert all(ψᵢ -> length(first(ψ)) == length(ψᵢ), ψ)

    # Output tensor
    ϕ = MPST(n)

    for j = 1:n
        ϕ[j] = ψ[1]
        for ψk in ψ[2:end]
            ϕ[j] += ψk[j]
        end
    end
    return ϕ
end
