function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:cutoff, :maxdim, :alg), <:Tuple{AbstractFloat, Integer, String}},typeof(product),MPO,MPS})   # time: 0.30468184
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{cutoff::Float64, maxdim::Int64, alg::String, nsweeps::Int64},typeof(product),MPO,MPS})   # time: 0.09260539
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{cutoff::Float64, maxdim::Int64, alg::String},typeof(product),MPO,MPS})   # time: 0.020713422
    isdefined(ITensors.ITensorMPS, Symbol("#94#109")) && Base.precompile(Tuple{getfield(ITensors.ITensorMPS, Symbol("#94#109")),Vector{Index{Int64}}})   # time: 0.003117791
    Base.precompile(Tuple{typeof(linkinds),MPO})   # time: 0.001600966
    Base.precompile(Tuple{typeof(linkdims),MPS})   # time: 0.001496124
end
