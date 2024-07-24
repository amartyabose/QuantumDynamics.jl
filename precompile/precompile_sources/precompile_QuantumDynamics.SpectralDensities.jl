function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(eval_spectrum),ExponentialCutoff,Real,Float64})   # time: 0.003293463
end
