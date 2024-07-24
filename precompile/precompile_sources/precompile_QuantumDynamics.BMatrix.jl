function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(get_B_matrix),Vector{Float64},Vector{Float64},Float64,Float64,Int64})   # time: 0.05446195
end
