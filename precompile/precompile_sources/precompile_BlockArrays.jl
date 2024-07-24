function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(getindex),Array{ComplexF64, 3},Block{1},Colon,Colon})   # time: 0.005674253
    Base.precompile(Tuple{typeof(getindex),Matrix{ComplexF64},Colon,Block{1}})   # time: 0.005599748
    Base.precompile(Tuple{typeof(getindex),Matrix{Float64},Block{1},Colon})   # time: 0.002814407
end
