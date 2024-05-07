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
    Base.precompile(Tuple{typeof(Transducers.AutoObjectsReStacker._getlength),Type{Tuple{Array{ComplexF64, 3}}}})
    Base.precompile(Tuple{typeof(Transducers.AutoObjectsReStacker._getlength),Type{Tuple{Float64}}})
    Base.precompile(Tuple{typeof(Transducers.AutoObjectsReStacker._getlength),Type{Tuple{Tuple{Float64}}}})
    Base.precompile(Tuple{typeof(Transducers.AutoObjectsReStacker._getlength),Type{Tuple{Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(transduce),IdentityTransducer,Function,InitOf{DefaultInitOf},Vector{Any},PreferParallel{@NamedTuple{simd::Val{false}}}})
    Base.precompile(Tuple{typeof(transduce),IdentityTransducer,Function,InitOf{DefaultInitOf},Vector{Vector{UInt64}},PreferParallel{@NamedTuple{simd::Val{false}}}})
    Base.precompile(Tuple{typeof(which(restack,(Any,)).generator.gen),Any,Any})
end
