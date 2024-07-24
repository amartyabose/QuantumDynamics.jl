function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(eigen),Matrix{ComplexF64}})   # time: 0.22581112
    Base.precompile(Tuple{typeof(generic_matmatmul!),Matrix{ComplexF64},Char,Char,Matrix{ComplexF64},Matrix{ComplexF64},MulAddMul{true, true, ComplexF64, ComplexF64}})   # time: 0.20332518
    Base.precompile(Tuple{typeof(inv),Matrix{ComplexF64}})   # time: 0.020965656
    Base.precompile(Tuple{typeof(*),Matrix{ComplexF64},Vector{ComplexF64},ComplexF64})   # time: 0.009077164
    Base.precompile(Tuple{typeof(diagm),Vector{Float64}})   # time: 0.00234337
end
