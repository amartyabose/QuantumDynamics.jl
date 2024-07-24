function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(copyat_or_push!),Vector{Vector{Matrix{ComplexF64}}},Int64,Vector{Matrix{ComplexF64}}})   # time: 0.002543469
end
