function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(_contract),DiagTensor{Float64, 2, Tuple{Index{Int64}, Index{Int64}}, Diag{Float64, Vector{Float64}}},DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.8315587
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.193965
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tensor{Number, 2, Combiner, Tuple{Index{Int64}, Index{Int64}}}})   # time: 0.1915788
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 6, NTuple{6, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.11647748
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 5, NTuple{5, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{Int64, 1, Tuple{Index{Int64}}, Dense{Int64, Vector{Int64}}}})   # time: 0.107398085
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.10537161
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 5, NTuple{5, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.10026999
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 5, NTuple{5, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.0917881
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.089302644
    Base.precompile(Tuple{Type{ITensor},Index{Int64},Index{Int64},Index{Int64},Any})   # time: 0.07781248
    Base.precompile(Tuple{typeof(_contract),DiagTensor{Float64, 2, Tuple{Index{Int64}, Index{Int64}}, Diag{Float64, Vector{Float64}}},DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.06912738
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.06854694
    Base.precompile(Tuple{Type{ITensor},Type{EmptyNumber},Index{Int64},Vararg{Index{Int64}}})   # time: 0.06462584
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{allow_alias::Bool},typeof(permute),ITensor,Index{Int64},Index{Int64}})   # time: 0.054129876
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tensor{Number, 3, Combiner, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})   # time: 0.04020583
    Base.precompile(Tuple{typeof(_permute),AllowAlias,DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tuple{Index{Int64}, Index{Int64}}})   # time: 0.03825284
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 6, NTuple{6, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.035376545
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.035146877
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.034502372
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.032824792
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{Int64, 1, Tuple{Index{Int64}}, Dense{Int64, Vector{Int64}}}})   # time: 0.032605045
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{tags::ITensors.TagSets.GenericTagSet{UInt256, 4}, positive::Bool},typeof(qr),ITensor,Vector{Index{Int64}},Vector{Index{Int64}}})   # time: 0.031952783
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.03127317
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Base.ReshapedArray{ComplexF64, 1, Adjoint{ComplexF64, Matrix{ComplexF64}}, Tuple{Base.MultiplicativeInverses.SignedMultiplicativeInverse{Int64}}}}},Tensor{Number, 2, Combiner, Tuple{Index{Int64}, Index{Int64}}}})   # time: 0.030868672
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 5, NTuple{5, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tensor{Number, 3, Combiner, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})   # time: 0.029218936
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{Float64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{Float64, Vector{Float64}}}})   # time: 0.024973504
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{Float64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{Float64, Vector{Float64}}}})   # time: 0.024660705
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{tags::ITensors.TagSets.GenericTagSet{UInt256, 4}, maxdim::Nothing},typeof(factorize),ITensor,Vector{Index{Int64}}})   # time: 0.022433827
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{Float64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{Float64, Vector{Float64}}}})   # time: 0.022244414
    Base.precompile(Tuple{typeof(combiner),Vector{Index{Int64}}})   # time: 0.021756263
    Base.precompile(Tuple{typeof(setindex!),ITensor,Number,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})   # time: 0.021435099
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tensor{Number, 3, Combiner, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})   # time: 0.021173552
    Base.precompile(Tuple{typeof(_contract),DenseTensor{Float64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{Float64, Vector{Float64}}},DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.020726543
    Base.precompile(Tuple{typeof(replaceind),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Index{Int64},Index{Int64}})   # time: 0.020467425
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tensor{Number, 2, Combiner, Tuple{Index{Int64}, Index{Int64}}}})   # time: 0.018494422
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tensor{Number, 3, Combiner, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})   # time: 0.015050119
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{mindim::Int64, maxdim::Int64, cutoff::Float64, ortho::String, which_decomp::Nothing, eigen_perturbation::Nothing, svd_alg::Nothing, tags::ITensors.TagSets.GenericTagSet{UInt256, 4}, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing, min_blockdim::Nothing},typeof(factorize),ITensor,Tuple{Index{Int64}, Index{Int64}}})   # time: 0.012996796
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{mindim::Int64, maxdim::Int64, cutoff::Float64, ortho::String, which_decomp::Nothing, eigen_perturbation::Nothing, svd_alg::Nothing, tags::ITensors.TagSets.GenericTagSet{UInt256, 4}, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing, min_blockdim::Nothing},typeof(factorize),ITensor,Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})   # time: 0.012296664
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{ishermitian::Bool, tags::ITensors.TagSets.GenericTagSet{UInt256, 4}, cutoff::Float64, maxdim::Int64, mindim::Int64},typeof(eigen),ITensor,Vector{Index{Int64}},Vector{Index{Int64}}})   # time: 0.011378994
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Index{Int64}},Vector{Index{Int64}}})   # time: 0.01028296
    Base.precompile(Tuple{typeof(prime),DenseTensor{Int64, 1, Tuple{Index{Int64}}, Dense{Int64, Vector{Int64}}}})   # time: 0.010067211
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Index{Int64}},Vector{Index{Int64}}})   # time: 0.008069372
    Base.precompile(Tuple{typeof(replaceinds),DiagTensor{Float64, 2, Tuple{Index{Int64}, Index{Int64}}, Diag{Float64, Vector{Float64}}},Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})   # time: 0.008042244
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Index{Int64}},Vector{Index{Int64}}})   # time: 0.007227665
    Base.precompile(Tuple{typeof(replaceinds),Tuple{Index{Int64}},Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})   # time: 0.006690711
    Base.precompile(Tuple{typeof(swapinds),ITensor,Vector{Index{Int64}},Vector{Index{Int64}}})   # time: 0.00654729
    Base.precompile(Tuple{typeof(setindex!),ITensor,Number,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})   # time: 0.006495713
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Base.ReshapedArray{ComplexF64, 1, Adjoint{ComplexF64, Matrix{ComplexF64}}, Tuple{Base.MultiplicativeInverses.SignedMultiplicativeInverse{Int64}}}}},Tensor{Number, 3, Combiner, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})   # time: 0.00635571
    Base.precompile(Tuple{typeof(setindex!),ITensor,Number,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, UInt64}})   # time: 0.005651247
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})   # time: 0.005479588
    Base.precompile(Tuple{typeof(swapinds),DenseTensor{Int64, 1, Tuple{Index{Int64}}, Dense{Int64, Vector{Int64}}},Vector{Index{Int64}},Vector{Index{Int64}}})   # time: 0.00545176
    Base.precompile(Tuple{typeof(replaceinds),ITensor,Vector{Pair{Index{Int64}, Index{Int64}}}})   # time: 0.005259279
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{dir::Nothing},typeof(combiner),Index{Int64}})   # time: 0.005232961
    Base.precompile(Tuple{typeof(replaceinds),NTuple{4, Index{Int64}},NTuple{4, Index{Int64}},NTuple{4, Index{Int64}}})   # time: 0.005178461
    Base.precompile(Tuple{typeof(replaceind),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Index{Int64},Index{Int64}})   # time: 0.005104084
    Base.precompile(Tuple{typeof(replaceinds),NTuple{4, Index{Int64}},Pair{Index{Int64}, Index{Int64}}})   # time: 0.005074502
    Base.precompile(Tuple{typeof(replaceinds),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},NTuple{4, Index{Int64}},NTuple{4, Index{Int64}}})   # time: 0.004824593
    Base.precompile(Tuple{Type{ITensor},Index{Int64},Index{Int64},Index{Int64},Index{Int64}})   # time: 0.004636257
    Base.precompile(Tuple{typeof(swapinds),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Index{Int64}},Vector{Index{Int64}}})   # time: 0.004480876
    Base.precompile(Tuple{typeof(indices),Index{Int64},Index{Int64},Vararg{Any}})   # time: 0.00419555
    Base.precompile(Tuple{typeof(swapinds),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Index{Int64}},Vector{Index{Int64}}})   # time: 0.004163788
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{dir::Arrow, tags::String},typeof(combiner),Index{Int64}})   # time: 0.004111328
    Base.precompile(Tuple{typeof(setindex!),ITensor,Number,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})   # time: 0.003795208
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{ortho::String, which_decomp::String},typeof(factorize),ITensor,Index{Int64}})   # time: 0.003759334
    Base.precompile(Tuple{typeof(setindex!),ITensor,Number,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, UInt64}})   # time: 0.003744295
    Base.precompile(Tuple{typeof(setindex!),ITensor,Int64,Pair{Index{Int64}, Int64}})   # time: 0.003428208
    Base.precompile(Tuple{typeof(itensor),DenseTensor{ComplexF64, 6, NTuple{6, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.003266663
    Base.precompile(Tuple{typeof(setindex!),ITensor,Int64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})   # time: 0.003047709
    Base.precompile(Tuple{typeof(replaceinds),ITensor,Pair{Index{Int64}, Index{Int64}}})   # time: 0.003012745
    Base.precompile(Tuple{typeof(setindex!),ITensor,AbstractArray,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})   # time: 0.002948997
    Base.precompile(Tuple{typeof(itensor),DenseTensor{ComplexF64, 5, NTuple{5, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})   # time: 0.002936247
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{mindim::Nothing, maxdim::Int64, cutoff::Float64, tags::ITensors.TagSets.GenericTagSet{UInt256, 4}, ortho::String, alg::Nothing, dir::Nothing, singular_values!::Nothing, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing, min_blockdim::Nothing},typeof(factorize_svd),ITensor,Index{Int64}})   # time: 0.00276233
    Base.precompile(Tuple{typeof(setindex!),ITensor,Float64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})   # time: 0.002727082
    Base.precompile(Tuple{typeof(dag),Vector{Index{Int64}}})   # time: 0.002404243
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{mindim::Nothing, maxdim::Int64, cutoff::Nothing, tags::ITensors.TagSets.GenericTagSet{UInt256, 4}, ortho::String, alg::Nothing, dir::Nothing, singular_values!::Nothing, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing, min_blockdim::Nothing},typeof(factorize_svd),ITensor,Index{Int64}})   # time: 0.002318838
    Base.precompile(Tuple{typeof(==),Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})   # time: 0.002318339
    Base.precompile(Tuple{typeof(setindex!),ITensor,AbstractArray,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})   # time: 0.002277125
    Base.precompile(Tuple{typeof(settags),DiagTensor{Float64, 2, Tuple{Index{Int64}, Index{Int64}}, Diag{Float64, Vector{Float64}}},ITensors.TagSets.GenericTagSet{UInt256, 4},Index{Int64}})   # time: 0.002247413
    Base.precompile(Tuple{typeof(setindex!),ITensor,ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})   # time: 0.002169709
    Base.precompile(Tuple{typeof(permute),Vector{Index{Int64}},Vector{Index{Int64}}})   # time: 0.002151957
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}},Vector{Index{Int64}}})   # time: 0.002060542
    Base.precompile(Tuple{typeof(combiner),Index{Int64},Index{Int64}})   # time: 0.00205604
    Base.precompile(Tuple{typeof(setindex!),ITensor,AbstractArray,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})   # time: 0.001838579
    Base.precompile(Tuple{Type{ITensor},Index{Int64},Index{Int64},Index{Int64}})   # time: 0.001819459
    Base.precompile(Tuple{typeof(settags),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},ITensors.TagSets.GenericTagSet{UInt256, 4},Index{Int64}})   # time: 0.00178896
    Base.precompile(Tuple{typeof(setindex!),ITensor,AbstractArray,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, UInt64}})   # time: 0.00155546
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{lefttags::ITensors.TagSets.GenericTagSet{UInt256, 4}, cutoff::Float64, maxdim::Int64},typeof(svd),ITensor,Vector{Index{Int64}}})   # time: 0.001504331
    Base.precompile(Tuple{typeof(setindex!),ITensor,AbstractArray,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, UInt64}})   # time: 0.001341752
    Base.precompile(Tuple{typeof(_setdiff),NTuple{4, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Vararg{Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})   # time: 0.001321041
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}},Vararg{Tuple{Index{Int64}, Index{Int64}}}})   # time: 0.001300464
    Base.precompile(Tuple{typeof(replaceinds),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Pair{Index{Int64}, Index{Int64}}})   # time: 0.001171704
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{mindim::Int64, maxdim::Int64, cutoff::Float64, tags::ITensors.TagSets.GenericTagSet{UInt256, 4}, ortho::String, eigen_perturbation::Nothing, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing},typeof(factorize_eigen),ITensor,Tuple{Index{Int64}, Index{Int64}}})   # time: 0.001154624
    Base.precompile(Tuple{typeof(settags),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},ITensors.TagSets.GenericTagSet{UInt256, 4},Index{Int64}})   # time: 0.001152628
    Base.precompile(Tuple{Type{Vector{IndexT} where IndexT<:Index},Index{Int64},Vararg{Index{Int64}}})   # time: 0.001142916
    Base.precompile(Tuple{typeof(replaceprime),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Pair{Int64, Int64}})   # time: 0.00111175
    Base.precompile(Tuple{Type{ITensor},Matrix{ComplexF64},Index{Int64},Index{Int64}})   # time: 0.00110208
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{Float64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{Float64, Vector{Float64}}},Vector{Pair{Index{Int64}, Index{Int64}}}})   # time: 0.001047001
    Base.precompile(Tuple{typeof(dag),AllowAlias,Tensor{Number, 3, Combiner, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})   # time: 0.001003416
    Base.precompile(Tuple{typeof(dag),AllowAlias,Tensor{Number, 2, Combiner, Tuple{Index{Int64}, Index{Int64}}}})   # time: 0.001001999
end
