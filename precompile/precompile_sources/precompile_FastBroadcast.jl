function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(fast_materialize),False,False,Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(muladd), Tuple{Float64, Matrix{ComplexF64}, Matrix{ComplexF64}}}})   # time: 0.003062036
    Base.precompile(Tuple{typeof(_rall),typeof(first),Tuple{Tuple{Bool, False}, Tuple{Bool, False}, Tuple{Bool, False}}})   # time: 0.001336542
end
