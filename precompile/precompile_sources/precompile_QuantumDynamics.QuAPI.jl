function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{Type{Path{ComplexF64}},Vector{Int64},Any,Int64})   # time: 0.002645703
    Base.precompile(Tuple{Type{QuAPIArgs}})   # time: 0.001081874
end
