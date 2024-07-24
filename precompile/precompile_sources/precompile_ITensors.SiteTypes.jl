function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{add_tags::String},typeof(siteinds),Int64,Int64})   # time: 0.008249542
    Base.precompile(Tuple{typeof(siteinds),Int64,Int64})   # time: 0.003714533
end
