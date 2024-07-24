function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    isdefined(Distributed, Symbol("#137#139")) && Base.precompile(Tuple{getfield(Distributed, Symbol("#137#139"))})   # time: 0.00208817
end
