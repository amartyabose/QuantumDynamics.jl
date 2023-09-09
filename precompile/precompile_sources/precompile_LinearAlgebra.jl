# Use
#    @warnpcfail precompile(args...)
# if you want to be warned when a precompile directive fails
macro warnpcfail(ex::Expr)
    modl = __module__
    file = __source__.file === nothing ? "?" : String(__source__.file)
    line = __source__.line
    quote
        $(esc(ex)) || @warn """precompile directive
     $($(Expr(:quote, ex)))
 failed. Please report an issue in $($modl) (after checking for duplicates) or remove this directive.""" _file=$file _line=$line
    end
end


function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(*),Matrix{ComplexF64},Vector{ComplexF64},ComplexF64})
    Base.precompile(Tuple{typeof(diagm),Vector{Float64}})
    Base.precompile(Tuple{typeof(gemm_wrapper!),Matrix{ComplexF64},Char,Char,Matrix{ComplexF64},Matrix{ComplexF64},MulAddMul{true, true, ComplexF64, ComplexF64}})
    Base.precompile(Tuple{typeof(mul!),Matrix{ComplexF64},Matrix{ComplexF64},Transpose{ComplexF64, Matrix{ComplexF64}},ComplexF64,ComplexF64})
    Base.precompile(Tuple{typeof(mul!),Matrix{ComplexF64},Transpose{ComplexF64, Matrix{ComplexF64}},Matrix{ComplexF64},ComplexF64,ComplexF64})
end
