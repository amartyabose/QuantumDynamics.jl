function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    isdefined(HDF5, Symbol("#134#137")) && Base.precompile(Tuple{getfield(HDF5, Symbol("#134#137"))})   # time: 0.001063876
    isdefined(HDF5, Symbol("#121#127")) && Base.precompile(Tuple{getfield(HDF5, Symbol("#121#127"))})   # time: 0.001052959
end
