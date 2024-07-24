module ComplexPISetup

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

function get_complex_time_propagator(Hamiltonian::AbstractMatrix{<:Complex}, β::Float64, t::Float64, N::Int64)
    Δtc = (t - 1im * β / 2) / N
    exp(-1im * Hamiltonian * Δtc), exp(1im * Hamiltonian * conj(Δtc))
end

end