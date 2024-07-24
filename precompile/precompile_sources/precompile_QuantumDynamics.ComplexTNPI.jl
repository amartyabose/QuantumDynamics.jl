function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(get_Bmat_MPO_left),Matrix{Float64},Vector,Int64,Vector{Matrix{ComplexF64}},Int64})   # time: 0.017964996
    Base.precompile(Tuple{typeof(get_Bmat_MPO_right),Matrix{Float64},Vector{Index{Int64}},Int64,Vector{Matrix{ComplexF64}},Int64})   # time: 0.015762536
    Base.precompile(Tuple{typeof(get_Bmat_MPO_left),Matrix{Float64},Vector{Index{Int64}},Int64,Vector{Matrix{ComplexF64}},Int64})   # time: 0.015465842
    Base.precompile(Tuple{typeof(get_Bmat_MPO_right),Matrix{Float64},Vector,Int64,Vector{Matrix{ComplexF64}},Int64})   # time: 0.008204346
end
