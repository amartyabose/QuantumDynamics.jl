module QuAPI

using Kronecker

using ..EtaCoefficients, ..SpectralDensities

struct States
    sbar :: Vector{Float64}
    Δs :: Vector{Float64}
    forward_ind :: Vector{UInt8}
    backward_ind :: Vector{UInt8}
end

struct Path{T}
    states :: Vector{UInt8}
    amplitude :: T
    sdim2 :: UInt8
end
function next_path(p::Path{T}) where {T<:Number}
    states = p.states
    sdim2 = p.sdim2
    for j=1:length(states)
        if states[j] < sdim2
            states[j] += 1
            break
        elseif states[j] == sdim2
            states[j] = 1
        end
    end
    return Path{T}(states, zero(T), sdim2), states != repeat([sdim2], length(states))
end

@inbounds function hash_path(states::Vector{UInt8}, sdim)
    factor = 1
    number = 0
    for s in states
        number += (s-1) * factor
        factor *= sdim
    end
    number + 1
end

@inbounds function unhash_path(path_num::Int, ntimes, sdim)
    path_num -= 1
    states = zeros(UInt8, ntimes+1)
    for j in 1:ntimes+1
        states[j] = path_num % sdim
        path_num = path_num ÷ sdim
    end
    states .+ 1
end

@inbounds function setup_simulation(H, ρ0, dt, η, svec, cutoff)
    sdim = size(H, 1)
    sdim_square = sdim^2
    U = exp(-1im * dt * H)
    fbU = U ⊗ conj(transpose(U))
    sbar = zeros(sdim_square)
    Δs = zeros(sdim_square)
    forward_ind = zeros(UInt8, sdim_square)
    backward_ind = zeros(UInt8, sdim_square)
    amplitudes = zeros(ComplexF64, sdim_square)
    count = 1
    for (forind, sf) in enumerate(svec)
        for (backind, sb) in enumerate(svec)
            sbar[count] = (sf + sb) / 2
            Δs[count] = sf - sb
            forward_ind[count] = forind
            backward_ind[count] = backind
            amplitudes[count] = ρ0[forind, backind]
            count += 1
        end
    end
    state_values = States(sbar, Δs, forward_ind, backward_ind)
    paths = Vector{Path{ComplexF64}}()
    for count = 1:sdim_square
        states = [count]
        if abs(amplitudes[count]) > cutoff
            push!(paths, Path{ComplexF64}(states, amplitudes[count] * exp(-state_values.Δs[count] * (real(η.η00) * state_values.Δs[count] + 2im * imag(η.η00) * state_values.sbar[count])), sdim_square))
        end
    end
    fbU, state_values, paths
end

@inbounds function get_influence(η::EtaCoefficients.EtaCoeffs, state_values::States, path::Vector{UInt8}, terminal::Bool, in_memory::Bool)
    Δsfinal = state_values.Δs[path[end]]
    if Δsfinal == 0
        return 1
    end
    num_time = length(path)
    if in_memory
        if terminal
            influence_exponent = real(η.η00) * Δsfinal + 2im * imag(η.η00) * state_values.sbar[path[end]]
            for i=2:num_time-1
                influence_exponent += real(η.η0m[num_time-i]) * state_values.Δs[path[i]] + 2im * imag(η.η0m[num_time-i]) * state_values.sbar[path[i]]
            end
            influence_exponent += real(η.η0e[num_time-1]) * state_values.Δs[path[1]] + 2im * imag(η.η0e[num_time-1]) * state_values.sbar[path[1]]
            return exp(-Δsfinal * influence_exponent)
        else
            influence_exponent = real(η.ηmm) * Δsfinal + 2im * imag(η.ηmm) * state_values.sbar[path[end]]
            for i=2:num_time-1
                influence_exponent += real(η.ηmn[num_time-i]) * state_values.Δs[path[i]] + 2im * imag(η.ηmn[num_time-i]) * state_values.sbar[path[i]]
            end
            influence_exponent += real(η.η0m[num_time-1]) * state_values.Δs[path[1]] + 2im * imag(η.η0m[num_time-1]) * state_values.sbar[path[1]]
            return exp(-Δsfinal * influence_exponent)
        end
    else
        if terminal
            influence_exponent = real(η.η00) * Δsfinal + 2im * imag(η.η00) * state_values.sbar[path[end]]
            for i=1:num_time-1
                influence_exponent += real(η.η0m[num_time-i]) * state_values.Δs[path[i]] + 2im * imag(η.η0m[num_time-i]) * state_values.sbar[path[i]]
            end
            return exp(-Δsfinal * influence_exponent)
        else
            influence_exponent = real(η.ηmm) * Δsfinal + 2im * imag(η.ηmm) * state_values.sbar[path[end]]
            for i=1:num_time-1
                influence_exponent += real(η.ηmn[num_time-i]) * state_values.Δs[path[i]] + 2im * imag(η.ηmn[num_time-i]) * state_values.sbar[path[i]]
            end
            return exp(-Δsfinal * influence_exponent)
        end
    end
end

@inbounds function propagate(;Hamiltonian, Jw::SpectralDensities.SpectralDensity, β::Real, ρ0, dt::Real, ntimes::Int, kmax::Int, cutoff=0.0, svec=[1.0, -1.0], verbose::Bool=false)
    if kmax > ntimes
        kmax = ntimes + 2
    end
    η = EtaCoefficients.calculate_η(Jw; β, dt, kmax)
    sdim = size(Hamiltonian, 1)
    ρs = zeros(ComplexF64, sdim, sdim, ntimes+1)
    ρs[:, :, 1] = ρ0
    fbU, state_values, paths = setup_simulation(Hamiltonian, ρ0, dt, η, svec, cutoff)

    sdim2 = sdim^2
    if verbose
        @info "Starting propagation within memory"
    end
    tmprho = zeros(ComplexF64, sdim^2)
    for i = 1:kmax
        if verbose
            @info "Step = $(i), #paths = $(length(paths))"
        end
        new_paths = Vector{Path{ComplexF64}}()
        for path in paths
            states = path.states
            tmprho .= zero(ComplexF64)
            tmprho[states[end]] = 1.0
            tmprho = fbU * tmprho * path.amplitude
            for (s, amp) in enumerate(tmprho)
                tmpstates = deepcopy(states)
                push!(tmpstates, s)
                ρs[state_values.forward_ind[s],state_values.backward_ind[s],i+1] += amp * get_influence(η, state_values, tmpstates, true, true)
                amplitude = amp * get_influence(η, state_values, tmpstates, false, true)
                if abs(amplitude) > cutoff
                    push!(new_paths, Path{ComplexF64}(tmpstates, amplitude, sdim2))
                end
            end
        end
        paths = new_paths
        new_paths = nothing
    end

    if verbose
        @info "Starting iteration"
    end
    path_max = repeat(Vector{UInt8}([sdim2]), length(paths[1].states)-1)
    max_length = hash_path(path_max, sdim2)
    checked = repeat([false], max_length)
    tmat = zeros(ComplexF64, max_length)
    for i = kmax+1:ntimes
        if verbose
            @info "Step = $(i), #paths = $(length(paths))"
        end
        tmat .= zero(ComplexF64)
        for path in paths
            trunc_path = path.states[2:end]
            tmat[hash_path(trunc_path, sdim2)] += path.amplitude
        end
        new_paths = Vector{Path{ComplexF64}}()
        checked .= false
        for path in paths
            states = path.states[2:end]
            hash_val = hash_path(states, sdim2)
            if checked[hash_val]
                continue
            else
                tmprho .= zero(ComplexF64)
                tmprho[states[end]] = 1.0
                tmprho = fbU * tmprho * tmat[hash_val]
                for (s, amp) in enumerate(tmprho)
                    tmpstates = deepcopy(states)
                    push!(tmpstates, s)
                    ρs[state_values.forward_ind[s],state_values.backward_ind[s],i+1] += amp * get_influence(η, state_values, tmpstates, true, true)
                    amplitude = amp * get_influence(η, state_values, tmpstates, false, true)
                    if abs(amplitude) > cutoff
                        push!(new_paths, Path{ComplexF64}(tmpstates, amplitude, sdim2))
                    end
                end
                checked[hash_val] = true
            end
        end
        paths = new_paths
        new_paths = nothing
    end
    ρs
end

@inbounds function get_path_amps(; η, state_values, path)
    interm_influence = zero(ComplexF64)
    i = length(path) - 1
    for sk = 2:i
        val = real(η.ηmm) * state_values.Δs[path[sk]] + 2im * imag(η.ηmm) * state_values.sbar[path[sk]]
        for skp = 2:sk-1
            val += real(η.ηmn[sk-skp]) * state_values.Δs[path[skp]] + 2im * imag(η.ηmn[sk-skp]) * state_values.sbar[path[skp]]
        end
        interm_influence += -state_values.Δs[path[sk]] * val
    end
    amplitude = exp(interm_influence)

    init_influence_0 = -state_values.Δs[path[1]] * (real(η.η00) * state_values.Δs[path[1]] + 2im * imag(η.η00) * state_values.sbar[path[1]])
    init_influence_m = -state_values.Δs[path[1]] * (real(η.ηmm) * state_values.Δs[path[1]] + 2im * imag(η.ηmm) * state_values.sbar[path[1]])
    for j = 2:i
        init_influence_0 += -state_values.Δs[path[j]] * (real(η.η0m[j-1] * state_values.Δs[path[1]]) + 2im * imag(η.η0m[j-1]) * state_values.sbar[path[1]])
        init_influence_m += -state_values.Δs[path[j]] * (real(η.ηmn[j-1] * state_values.Δs[path[1]]) + 2im * imag(η.ηmn[j-1]) * state_values.sbar[path[1]])
    end
    init_influence_0 = exp(init_influence_0)
    init_influence_m = exp(init_influence_m)

    final_influence_0 = -state_values.Δs[path[end]] * (real(η.η00) * state_values.Δs[path[end]] + 2im * imag(η.η00) * state_values.sbar[path[end]])
    final_influence_m = -state_values.Δs[path[end]] * (real(η.ηmm) * state_values.Δs[path[end]] + 2im * imag(η.ηmm) * state_values.sbar[path[end]])
    for j = 2:i
        final_influence_0 += -state_values.Δs[path[end]] * (real(η.η0m[i+1-j] * state_values.Δs[path[j]]) + 2im * imag(η.η0m[i+1-j]) * state_values.sbar[path[j]])
        final_influence_m += -state_values.Δs[path[end]] * (real(η.ηmn[i+1-j] * state_values.Δs[path[j]]) + 2im * imag(η.ηmn[i+1-j]) * state_values.sbar[path[j]])
    end
    final_influence_0 = exp(final_influence_0)
    final_influence_m = exp(final_influence_m)

    term_influence_0e = exp(-state_values.Δs[path[end]] * (real(η.η0e[i]) * state_values.Δs[path[1]] + 2im * imag(η.η0e[i]) * state_values.sbar[path[1]]))
    term_influence_me = exp(-state_values.Δs[path[end]] * (real(η.η0m[i]) * state_values.Δs[path[1]] + 2im * imag(η.η0m[i]) * state_values.sbar[path[1]]))
    term_influence_0m = exp(-state_values.Δs[path[end]] * (real(η.η0m[i]) * state_values.Δs[path[1]] + 2im * imag(η.η0m[i]) * state_values.sbar[path[1]]))
    term_influence_mn = exp(-state_values.Δs[path[end]] * (real(η.ηmn[i]) * state_values.Δs[path[1]] + 2im * imag(η.ηmn[i]) * state_values.sbar[path[1]]))

    amplitude * init_influence_0 * final_influence_0 * term_influence_0e, amplitude * init_influence_0 * final_influence_m * term_influence_0m, amplitude * init_influence_m * final_influence_0 * term_influence_me, amplitude * init_influence_m * final_influence_m * term_influence_mn
end

@inbounds function build_propagator(;Hamiltonian, Jw::SpectralDensities.SpectralDensity, β::Real, dt::Real, ntimes::Int, cutoff=0.0, svec=[1.0, -1.0], verbose::Bool=false)
    η = EtaCoefficients.calculate_η(Jw; β, dt, kmax=ntimes)
    sdim = size(Hamiltonian, 1)
    sdim2 = sdim^2
    fbU, state_values, _ = setup_simulation(Hamiltonian, ones(sdim, sdim), dt, η, svec, cutoff)

    if verbose
        @info "Starting propagation within memory"
    end
    U0e = zeros(ComplexF64, sdim2, sdim2, ntimes)
    U0m = zeros(ComplexF64, sdim2, sdim2, ntimes)
    Ume = zeros(ComplexF64, sdim2, sdim2, ntimes)
    Umn = zeros(ComplexF64, sdim2, sdim2, ntimes)
    for i = 1:ntimes
        if verbose
            @info "Step = $(i)"
        end
        for path_num = 1:sdim2^(i+1)
            states = unhash_path(path_num, i, sdim2)
            state_pairs = collect(zip(states, states[2:end]))
            bare_amplitude = prod([fbU[s[1], s[2]] for s in state_pairs])
            if abs(bare_amplitude) < cutoff
                continue
            end
            amplitudes = get_path_amps(; η, state_values, path=states)
            U0e[states[end], states[1], i] += amplitudes[1] * bare_amplitude
            U0m[states[end], states[1], i] += amplitudes[2] * bare_amplitude
            Ume[states[end], states[1], i] += amplitudes[3] * bare_amplitude
            Umn[states[end], states[1], i] += amplitudes[4] * bare_amplitude
        end
    end
    U0e, U0m, Ume, Umn
end

end