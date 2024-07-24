function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{Colon,Int64,Static.StaticInt{U}})   # time: 0.001697126
end
