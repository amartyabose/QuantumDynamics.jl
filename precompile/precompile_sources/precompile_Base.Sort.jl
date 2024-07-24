function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{rev::Bool, by::typeof(abs)},typeof(sortperm),Vector{Float64}})   # time: 0.015998332
end
