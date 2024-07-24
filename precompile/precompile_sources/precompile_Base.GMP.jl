function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(prod),Vector{BigInt}})   # time: 0.003046785
end
