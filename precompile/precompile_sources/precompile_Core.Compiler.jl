function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(apply_type_tfunc),InferenceLattice{PartialsLattice{ConstsLattice}},Any,Any,Vararg{Any}})   # time: 0.002964995
end
