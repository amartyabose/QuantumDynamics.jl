function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(pairs),@NamedTuple{simd::Base.Val{false}}})   # time: 0.001312626
end
