function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(get_complex_time_propagator),Matrix{ComplexF64},Float64,Float64,Int64})   # time: 0.006184968
end
