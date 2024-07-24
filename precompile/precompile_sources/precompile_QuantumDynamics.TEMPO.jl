function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    isdefined(QuantumDynamics.TEMPO, Symbol("#2#9")) && Base.precompile(Tuple{getfield(QuantumDynamics.TEMPO, Symbol("#2#9")),Int64})   # time: 0.00143521
end
