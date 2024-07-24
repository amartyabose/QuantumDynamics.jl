function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(svd_recursive),Matrix{ComplexF64}})   # time: 0.3206468
    Base.precompile(Tuple{typeof(permutedims!),Exposed{Array{ComplexF64, 3}, Array{ComplexF64, 3}},Exposed{Array{ComplexF64, 3}, Array{ComplexF64, 3}},Tuple{Int64, Int64, Int64}})   # time: 0.16816528
    Base.precompile(Tuple{typeof(_contract!),Matrix{ComplexF64},Matrix{ComplexF64},Matrix{ComplexF64},ContractionProperties{2, 2, 2},Int64,Int64})   # time: 0.1125149
    Base.precompile(Tuple{typeof(permutedims!),Exposed{Array{ComplexF64, 5}, Array{ComplexF64, 5}},Exposed{Array{ComplexF64, 5}, Array{ComplexF64, 5}},NTuple{5, Int64}})   # time: 0.09195845
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{mindim::Int64, maxdim::Int64, cutoff::Float64, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing},typeof(truncate!!),AbstractArray})   # time: 0.08258921
    Base.precompile(Tuple{typeof(permutedims!),Exposed{Array{ComplexF64, 4}, Array{ComplexF64, 4}},Exposed{Array{ComplexF64, 4}, Array{ComplexF64, 4}},NTuple{4, Int64}})   # time: 0.06814345
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 5},Array{ComplexF64, 3},Array{ComplexF64, 4},ContractionProperties{3, 4, 5},Int64,Int64})   # time: 0.054376498
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 6},Array{ComplexF64, 5},Array{ComplexF64, 3},ContractionProperties{5, 3, 6},Int64,Int64})   # time: 0.05248866
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 4},Array{ComplexF64, 5},Vector{ComplexF64},ContractionProperties{5, 1, 4},Int64,Int64})   # time: 0.044378538
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 3},Matrix{ComplexF64},Array{ComplexF64, 3},ContractionProperties{2, 3, 3},Int64,Int64})   # time: 0.03499588
    Base.precompile(Tuple{typeof(mul!!),Matrix{ComplexF64},Transpose{ComplexF64, Matrix{ComplexF64}},Matrix{ComplexF64},ComplexF64,ComplexF64})   # time: 0.025217712
    Base.precompile(Tuple{typeof(similartype),Type{Vector{ComplexF64}},Tuple{Int64, Int64}})   # time: 0.023963083
    Base.precompile(Tuple{typeof(_contract_scalar_perm!),ReshapedArray{ComplexF64, _A, Vector{ComplexF64}, Tuple{}} where _A,ReshapedArray{ComplexF64, _A, Vector{ComplexF64}, Tuple{}} where _A,Any,Float64,Int64})   # time: 0.01690239
    Base.precompile(Tuple{typeof(mul!!),Matrix{ComplexF64},Matrix{ComplexF64},Transpose{ComplexF64, Matrix{ComplexF64}},ComplexF64,ComplexF64})   # time: 0.01507414
    Base.precompile(Tuple{typeof(mul!!),Matrix{ComplexF64},Transpose{ComplexF64, Matrix{ComplexF64}},Transpose{ComplexF64, Matrix{ComplexF64}},ComplexF64,ComplexF64})   # time: 0.014442757
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{3, 4, 5}})   # time: 0.014367547
    Base.precompile(Tuple{typeof(cpu),Vector{Float64}})   # time: 0.0120354
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{6}},NTuple{5, Int64},Tuple{Int64, Int64, Int64}})   # time: 0.006619708
    Base.precompile(Tuple{typeof(_contract_scalar_perm!),ReshapedArray{ComplexF64, _A, Vector{ComplexF64}, Tuple{}} where _A,ReshapedArray{Float64, _A, Vector{Float64}, Tuple{}} where _A,Any,ComplexF64,Int64})   # time: 0.00633454
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{5}},Tuple{Int64, Int64, Int64},NTuple{4, Int64}})   # time: 0.00552642
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{alg::LinearAlgebra.QRIteration},typeof(svd_catch_error),Matrix{ComplexF64}})   # time: 0.004875454
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{2}},Tuple{Int64, Int64},Tuple{Int64, Int64}})   # time: 0.004044584
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{3}},Tuple{Int64, Int64},Tuple{Int64, Int64, Int64}})   # time: 0.003744539
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{mindim::Nothing, maxdim::Int64, cutoff::Nothing, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing},typeof(truncate!!),Vector{Float64}})   # time: 0.003592702
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{alg::LinearAlgebra.DivideAndConquer},typeof(svd_catch_error),Matrix{ComplexF64}})   # time: 0.003041628
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{mindim::Nothing, maxdim::Int64, cutoff::Float64, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing},typeof(truncate!!),Vector{Float64}})   # time: 0.002718122
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{4, 3, 3}})   # time: 0.002474662
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{6, 3, 5}})   # time: 0.002456166
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{5, 3, 4}})   # time: 0.002443657
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{6, 4, 6}})   # time: 0.00241938
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 4},Array{ComplexF64, 4},Array{ComplexF64, 4},ContractionProperties{4, 4, 4},Int64,Int64})   # time: 0.002407164
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{5, 3, 6}})   # time: 0.002395832
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 4},Array{ComplexF64, 4},Matrix{ComplexF64},ContractionProperties{4, 2, 4},Int64,Int64})   # time: 0.002343499
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{2, 3, 3}})   # time: 0.002343462
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{3, 2, 3}})   # time: 0.002326959
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{3, 4, 3}})   # time: 0.002326793
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{4, 2, 4}})   # time: 0.002325956
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 4},Array{ComplexF64, 3},Array{ComplexF64, 3},ContractionProperties{3, 3, 4},Int64,Int64})   # time: 0.002323542
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{4, 3, 5}})   # time: 0.002321127
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{5, 4, 5}})   # time: 0.002309367
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 6},Array{ComplexF64, 6},Array{ComplexF64, 4},ContractionProperties{6, 4, 6},Int64,Int64})   # time: 0.002283915
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 4},Array{ComplexF64, 4},Array{ComplexF64, 4},ContractionProperties{4, 4, 4},Bool,Bool})   # time: 0.002265375
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{5, 1, 4}})   # time: 0.002241915
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 4},Array{ComplexF64, 5},Array{ComplexF64, 3},ContractionProperties{5, 3, 4},Int64,Int64})   # time: 0.002234334
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{4, 1, 3}})   # time: 0.002224918
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 5},Array{ComplexF64, 4},Array{ComplexF64, 3},ContractionProperties{4, 3, 5},Int64,Int64})   # time: 0.002217999
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 5},Array{ComplexF64, 6},Array{ComplexF64, 3},ContractionProperties{6, 3, 5},Int64,Int64})   # time: 0.002188623
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 6},Array{ComplexF64, 6},Array{ComplexF64, 4},ContractionProperties{6, 4, 6},Bool,Bool})   # time: 0.002175664
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 4},Array{ComplexF64, 3},Array{ComplexF64, 3},ContractionProperties{3, 3, 4},Bool,Bool})   # time: 0.002157623
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{3, 3, 4}})   # time: 0.002156583
    Base.precompile(Tuple{typeof(*),Bool,Dense{ComplexF64, Vector{ComplexF64}}})   # time: 0.002155583
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 4},Array{ComplexF64, 4},Matrix{ComplexF64},ContractionProperties{4, 2, 4},Bool,Bool})   # time: 0.002130956
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 5},Array{ComplexF64, 5},Array{ComplexF64, 4},ContractionProperties{5, 4, 5},Int64,Int64})   # time: 0.002108751
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 3},Array{ComplexF64, 4},Vector{ComplexF64},ContractionProperties{4, 1, 3},Int64,Int64})   # time: 0.002105251
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 4},Array{ComplexF64, 5},Array{ComplexF64, 3},ContractionProperties{5, 3, 4},Bool,Bool})   # time: 0.002081707
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 6},Array{ComplexF64, 5},Array{ComplexF64, 3},ContractionProperties{5, 3, 6},Bool,Bool})   # time: 0.002080167
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 5},Array{ComplexF64, 3},Array{ComplexF64, 4},ContractionProperties{3, 4, 5},Bool,Bool})   # time: 0.002071916
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{4}},NTuple{5, Int64},Tuple{Int64}})   # time: 0.002065665
    Base.precompile(Tuple{typeof(_contract_scalar_perm!),ReshapedArray{ComplexF64, _A, Vector{ComplexF64}, Tuple{}} where _A,ReshapedArray{ComplexF64, _A, Vector{ComplexF64}, Tuple{}} where _A,Any,Float64,Bool})   # time: 0.002061923
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 3},Array{ComplexF64, 4},Array{ComplexF64, 3},ContractionProperties{4, 3, 3},Bool,Bool})   # time: 0.002050292
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 5},Array{ComplexF64, 4},Array{ComplexF64, 3},ContractionProperties{4, 3, 5},Bool,Bool})   # time: 0.002049918
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{3, 3, 2}})   # time: 0.00204258
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 5},Array{ComplexF64, 6},Array{ComplexF64, 3},ContractionProperties{6, 3, 5},Bool,Bool})   # time: 0.002040793
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 3},Array{ComplexF64, 4},Array{ComplexF64, 3},ContractionProperties{4, 3, 3},Int64,Int64})   # time: 0.002002792
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 3},Array{ComplexF64, 3},Matrix{ComplexF64},ContractionProperties{3, 2, 3},Bool,Bool})   # time: 0.001992166
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 5},Array{ComplexF64, 5},Array{ComplexF64, 4},ContractionProperties{5, 4, 5},Bool,Bool})   # time: 0.001982043
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 3},Array{ComplexF64, 3},Array{ComplexF64, 4},ContractionProperties{3, 4, 3},Int64,Int64})   # time: 0.001968292
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 3},Array{ComplexF64, 4},Vector{ComplexF64},ContractionProperties{4, 1, 3},Bool,Bool})   # time: 0.00196754
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 4},Array{ComplexF64, 5},Vector{ComplexF64},ContractionProperties{5, 1, 4},Bool,Bool})   # time: 0.001966043
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 3},Array{ComplexF64, 3},Matrix{ComplexF64},ContractionProperties{3, 2, 3},Int64,Int64})   # time: 0.001946457
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 3},Matrix{ComplexF64},Array{ComplexF64, 3},ContractionProperties{2, 3, 3},Bool,Bool})   # time: 0.001924669
    Base.precompile(Tuple{typeof(_contract!),Matrix{ComplexF64},Array{ComplexF64, 3},Array{ComplexF64, 3},ContractionProperties{3, 3, 2},Int64,Int64})   # time: 0.001841252
    Base.precompile(Tuple{typeof(_contract!),Array{ComplexF64, 3},Array{ComplexF64, 3},Array{ComplexF64, 4},ContractionProperties{3, 4, 3},Bool,Bool})   # time: 0.001831331
    Base.precompile(Tuple{typeof(_contract!),Matrix{ComplexF64},Matrix{ComplexF64},Matrix{ComplexF64},ContractionProperties{2, 2, 2},Bool,Bool})   # time: 0.001786622
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{4, 4, 4}})   # time: 0.001711414
    Base.precompile(Tuple{typeof(_contract!),Matrix{ComplexF64},Array{ComplexF64, 3},Array{ComplexF64, 3},ContractionProperties{3, 3, 2},Bool,Bool})   # time: 0.001612835
    Base.precompile(Tuple{Type{ContractionProperties},Tuple{Int64, Int64},Tuple{Int64, Int64},Tuple{Int64, Int64}})   # time: 0.001580334
    Base.precompile(Tuple{typeof(compute_perms!),ContractionProperties{2, 2, 2}})   # time: 0.001516163
end
