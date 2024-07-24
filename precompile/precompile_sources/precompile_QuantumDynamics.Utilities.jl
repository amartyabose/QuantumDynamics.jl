function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(commutator),Matrix{ComplexF64},Matrix{Float64}})   # time: 0.09811818
    Base.precompile(Tuple{typeof(build_path_amplitude_mps),Matrix{ComplexF64},Vector{Index{Int64}}})   # time: 0.039044563
    Base.precompile(Tuple{typeof(apply_contract_propagator),MPS,MPO})   # time: 0.027084377
    Base.precompile(Tuple{typeof(convert_ITensor_to_matrix),Any,Index{Int64},Index{Int64}})   # time: 0.012738292
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{discrete::Bool},typeof(trapezoid_alg),Execution{:seq},Vector{Float64},Vector{ComplexF64}})   # time: 0.009124904
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{discrete::Bool},typeof(trapezoid_alg),Execution{:seq},Vector{Float64},Vector{Float64}})   # time: 0.007852462
    Base.precompile(Tuple{typeof(unhash_path),Int64,Int64,Int64})   # time: 0.007355916
    Base.precompile(Tuple{typeof(path_amplitude_to_propagator),MPS})   # time: 0.006452293
    Base.precompile(Tuple{typeof(extend_path_amplitude_mps_beyond_memory),MPS,Matrix{ComplexF64},Vector{Index{Int64}}})   # time: 0.004301375
    Base.precompile(Tuple{typeof(has_small_changes),Any,Int64})   # time: 0.001087
end
