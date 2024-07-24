function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(hcat),Float64,Float64})   # time: 0.005790917
    Base.precompile(Tuple{typeof(hvcat),Tuple{Int64, Int64},ComplexF64,Float64,Vararg{Float64}})   # time: 0.001430542
end
