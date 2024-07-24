function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(func_BRME),Matrix{ComplexF64},Params,Float64})   # time: 0.062367275
    Base.precompile(Tuple{typeof(func_BRME),Any,Params,Float64})   # time: 0.035462204
    Base.precompile(Tuple{typeof(func_BRME),Matrix{ComplexF64},Any,Any})   # time: 0.003775296
end
