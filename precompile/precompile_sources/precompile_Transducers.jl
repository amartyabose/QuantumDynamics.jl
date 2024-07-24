function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(_transduce_assoc_nocomplete),AbstractReduction,InitOf{DefaultInitOf},Any,Any,NoopDACContext})   # time: 0.038440675
    Base.precompile(Tuple{typeof(transduce),IdentityTransducer,Function,InitOf{DefaultInitOf},Vector{Vector{UInt64}},ThreadedEx{@NamedTuple{simd::Val{false}}}})   # time: 0.021618038
    Base.precompile(Tuple{typeof(transduce),IdentityTransducer,Function,InitOf{DefaultInitOf},Vector{Any},ThreadedEx{@NamedTuple{simd::Val{false}}}})   # time: 0.020411786
    Base.precompile(Tuple{typeof(maybe_set_simd),ThreadedEx{@NamedTuple{}},Val{false}})   # time: 0.006197332
    Base.precompile(Tuple{typeof(maybe_set_simd),SequentialEx{@NamedTuple{}},Val{false}})   # time: 0.004606373
    Base.precompile(Tuple{typeof(__reduce_dummy),Function,InitOf{DefaultInitOf},SizedReducible{SubArray{Any, 1, Vector{Any}, Tuple{UnitRange{Int64}}, true}, Int64}})   # time: 0.001086083
    Base.precompile(Tuple{typeof(__reduce_dummy),Function,InitOf{DefaultInitOf},SizedReducible{SubArray{Vector{UInt64}, 1, Vector{Vector{UInt64}}, Tuple{UnitRange{Int64}}, true}, Int64}})   # time: 0.00107679
end
