function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{stale_age::Int64},typeof(trymkpidlock),Function,Vararg{Any}})   # time: 0.016402388
end
