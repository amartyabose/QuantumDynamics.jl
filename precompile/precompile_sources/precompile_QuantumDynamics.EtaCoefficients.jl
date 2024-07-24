function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(calculate_η),Vector{Float64},Vector{Float64},Float64,Float64,Int64,Bool,Bool,Bool})   # time: 0.42108637
end
