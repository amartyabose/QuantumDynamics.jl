function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(muladd),Float64,Matrix{ComplexF64},Matrix{ComplexF64}})   # time: 0.016755704
end
