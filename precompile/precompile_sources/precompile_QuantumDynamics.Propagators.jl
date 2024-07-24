function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{Hamiltonian::Matrix{ComplexF64}, dt::Float64, ntimes::Int64},typeof(calculate_bare_propagators)})   # time: 0.4211968
end
