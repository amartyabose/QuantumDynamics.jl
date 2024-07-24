function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(calculate_residuals),Matrix{ComplexF64},Matrix{ComplexF64},Matrix{ComplexF64},Float64,Float64,Function,Float64})   # time: 0.002866798
    Base.precompile(Tuple{typeof(ODE_DEFAULT_NORM),Matrix{ComplexF64},Float64})   # time: 0.002536167
end
