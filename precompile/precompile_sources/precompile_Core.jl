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
    Base.precompile(Tuple{Type{Array{Symbol}},UndefInitializer,Int64})
    Base.precompile(Tuple{Type{NamedTuple{(:algorithm,)}},Tuple{String}})
    Base.precompile(Tuple{Type{NamedTuple{(:ξ, :ωc)}},Tuple{Float64, Float64}})
    Base.precompile(Tuple{Type{NamedTuple{(:ϵ, :Δ)}},Tuple{Float64, Float64}})
    Base.precompile(Tuple{Type{Union{}},Float64})
    Base.precompile(Tuple{typeof(Core.Compiler.eltype),Core.Type{Base.Vector{Base.Vector{Base.Vector{Core.UInt64}}}}})
    Base.precompile(Tuple{typeof(Core.Compiler.eltype),Core.Type{Base.Vector{Base.Vector{Core.Any}}}})
    Base.precompile(Tuple{typeof(Core.Compiler.eltype),Core.Type{Base.Vector{Base.Vector}}})
    Base.precompile(Tuple{typeof(Core.Compiler.eltype),Core.Type{Base.Vector{Tuple{Core.Int64, Core.Int64}}}})
    Base.precompile(Tuple{typeof(Core.Compiler.setindex!),Base.Vector{Core.Compiler.CallInfo},Core.Compiler.InvokeCallInfo,Int64})
end
