function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(ode_interpolant),Float64,Float64,Matrix{ComplexF64},Matrix{ComplexF64},Vector{Matrix{ComplexF64}},Tsit5ConstantCache,Nothing,Type{Val{0}},Nothing})   # time: 0.003963668
    Base.precompile(Tuple{typeof(default_controller),Tsit5{typeof(trivial_limiter!), typeof(trivial_limiter!), False},Tsit5ConstantCache,Rational{Int64},Nothing,Nothing})   # time: 0.002928082
    Base.precompile(Tuple{typeof(gamma_default),Tsit5{typeof(trivial_limiter!), typeof(trivial_limiter!), False}})   # time: 0.002402417
end
