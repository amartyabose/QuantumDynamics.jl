using LinearAlgebra
using FLoops

"""
    ExecType{exec}
This empty struct holds the execution paradigm of various algorithms. `exec` currently can be one of `seq` and `par`.
"""
struct Execution{calculation} end
Execution(c) = Execution{Symbol(c)}
macro Execution_str(c)
    return :(Execution{$(Expr(:quote, Symbol(c)))})
end

"""
    get_BLAS_implementation()
Reports the BLAS implementation under use. The default implementation used by Julia is OpenBlas. MKL is supported through the external package MKL.jl, which needs to be installed and loaded before the loading of QuantumDynamics.jl
"""
get_BLAS_implementation() = BLAS.get_config()

trapezoid_alg(::Execution"seq", x, y; discrete::Bool=false) = discrete ? sum(y) : sum((y[1:end-1] .+ y[2:end]) .* (x[2:end] .- x[1:end-1])) / 2
function trapezoid_alg(::Execution"par", x, y; discrete::Bool=false)
    if discrete
        @floop for val in y
            @reduce ans = zero(eltype(y)) + val
        end
        ans
    else
        len = length(y)
        @floop for j = 1:len-1
            @reduce ans = zero(eltype(y)) + (x[j+1] - x[j]) * (y[j+1] + y[j]) / 2
        end
        ans
    end
end

"""
    trapezoid(x, y; discrete::Bool=false, exec="seq")
Returns the trapezoidal integration of y with respect to x. If discrete is set to `true`, then returns sum of y. `exec` encodes the execution paradigm and is one of `seq` or `par`.
"""
trapezoid(x, y; discrete::Bool=false, exec="seq") = trapezoid_alg(Execution(exec)(), x, y; discrete)

@doc raw"""
    fourier_transform(time::AbstractArray{<:Real}, corr::AbstractArray{<:Complex}; full=true, unitary=false)
Returns the Fourier transform of `corr`:

``C(\omega) = \mathcal{N}\,\int_a^\infty C(t)\,e^{i\omega t}\,dt``

If `unitary` is true, then ``\mathcal{N}=\frac{1}{2\pi}``, else ``\mathcal{N}=1``.

If `full` is true, then ``a=-\infty``, else ``a=0``.
"""
function fourier_transform(time::AbstractArray{<:Real}, corr::AbstractArray{<:Complex}; full=true, unitary=false, ωmin::Union{Nothing,Float64}=nothing, ωmax::Union{Nothing,Float64}=nothing, flip_sign=false)
    dt = time[2] - time[1]
    ωlim = π / dt
    dω = π / time[end]
    ωm = isnothing(ωmin) ? -ωlim : max(-ωlim, ωmin)
    ωM = isnothing(ωmax) ? ωlim : min(ωlim, ωmax)
    ω = ωm:dω:ωM |> collect
    spect = zeros(Complex{real(eltype(corr))}, length(ω))
    sign = flip_sign ? -1im : 1im
    if full
        for (l, w) in enumerate(ω)
            spect[l] = Utilities.trapezoid(time, corr .* exp.(sign * w * time) + conj.(corr) .* exp.(-sign * w * time))
        end
    else
        for (l, w) in enumerate(ω)
            spect[l] = Utilities.trapezoid(time, corr .* exp.(sign * w * time))
        end
    end

    spect ./= unitary ? sqrt(2π) : 1

    ω, spect
end

function inverse_fourier_transform(ω::AbstractArray{<:Real}, spect::AbstractArray{<:Complex}; unitary=false)
    dω = ω[2] - ω[1]
    dt = π / ω[end]
    tmax = π / dω
    time = 0:dt:tmax
    data = zeros(Complex{real(eltype(spect))}, length(time))
    for (l, t) in enumerate(time)
        data[l] = Utilities.trapezoid(ω, spect .* exp.(-1im * ω * t))
    end

    data ./= unitary ? sqrt(2π) : 2π
    time, data
end


function finite_difference_coeffs(n, order)
    M = zeros(order + 1, order + 1)
    M[1, :] .= 1
    for j = 2:order+1
        for k = 2:order+1
            M[j, k] = (-k + 1)^(j - 1)
        end
    end

    v = zeros(order + 1)
    v[n+1] = 1

    M \ v
end