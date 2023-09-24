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


const __bodyfunction__ = Dict{Method,Any}()

# Find keyword "body functions" (the function that contains the body
# as written by the developer, called after all missing keyword-arguments
# have been assigned values), in a manner that doesn't depend on
# gensymmed names.
# `mnokw` is the method that gets called when you invoke it without
# supplying any keywords.
function __lookup_kwbody__(mnokw::Method)
    function getsym(ast, arg)
        isa(arg, Symbol) && return arg
        isa(arg, GlobalRef) && return arg.name
        if isa(arg, Core.SSAValue)
            arg = ast.code[arg.id]
            return getsym(ast, arg)
        end
    end

    f = get(__bodyfunction__, mnokw, nothing)
    if f === nothing
        fmod = mnokw.module
        # The lowered code for `mnokw` should look like
        #   %1 = mkw(kwvalues..., #self#, args...)
        #        return %1
        # where `mkw` is the name of the "active" keyword body-function.
        ast = Base.uncompressed_ast(mnokw)
        if isa(ast, Core.CodeInfo) && length(ast.code) >= 2
            callexpr = ast.code[end-1]
            if isa(callexpr, Expr) && callexpr.head == :call
                fsym = callexpr.args[1]
                if isa(fsym, Symbol)
                    f = getfield(fmod, fsym)
                elseif isa(fsym, GlobalRef)
                    if fsym.mod === Core && fsym.name === :_apply
                        f = getfield(mnokw.module, getsym(ast, callexpr.args[2]))
                    elseif fsym.mod === Core && fsym.name === :_apply_iterate
                        f = getfield(mnokw.module, getsym(ast, callexpr.args[3]))
                    else
                        f = getfield(fsym.mod, fsym.name)
                    end
                else
                    f = missing
                end
            else
                f = missing
            end
        else
            f = missing
        end
        __bodyfunction__[mnokw] = f
    end
    return f
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{Type{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}}},typeof(real),Tuple{Vector{Float64}}})
    Base.precompile(Tuple{Type{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}}},typeof(*),Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(/), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(-), Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}}, Matrix{Float64}}}, Float64}})
    Base.precompile(Tuple{Type{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}}},typeof(*),Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(/), Tuple{Matrix{ComplexF64}, Matrix{Float64}}}, Float64}})
    Base.precompile(Tuple{Type{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}}},typeof(/),Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(-), Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}}, Matrix{Float64}}})
    Base.precompile(Tuple{Type{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}}},typeof(muladd),Tuple{Float64, Matrix{ComplexF64}, Matrix{ComplexF64}}})
    Base.precompile(Tuple{Type{Fix1},Type{MappingRF},Type})
    Base.precompile(Tuple{Type{Matrix{ComplexF64}},Matrix{Float64}})
    Base.precompile(Tuple{typeof(==),Matrix{ComplexF64},Matrix{ComplexF64}})
    Base.precompile(Tuple{typeof(Base.Cartesian._nloops),Int64,Symbol,Symbol,Expr,Vararg{Expr}})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:init,), Tuple{Int64}},typeof(mapreduce),Type,Function,Tuple{Int64}})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:length,), Tuple{Int64}},typeof(range),Float64,Float64})
    Base.precompile(Tuple{typeof(_getindex),IndexLinear,Array{ComplexF64, 3},Int64,Slice{OneTo{Int64}},Vararg{Slice{OneTo{Int64}}}})
    Base.precompile(Tuple{typeof(any),Function,Matrix{ComplexF64}})
    Base.precompile(Tuple{typeof(broadcasted),Function,Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(-), Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}},Matrix{Float64}})
    Base.precompile(Tuple{typeof(broadcasted),Function,Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(/), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(-), Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}}, Matrix{Float64}}},Float64})
    Base.precompile(Tuple{typeof(broadcasted),Function,Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(/), Tuple{Matrix{ComplexF64}, Matrix{Float64}}},Float64})
    Base.precompile(Tuple{typeof(broadcasted),Function,Base.Broadcast.Broadcasted{Base.Broadcast.Style{Tuple}, Nothing, typeof(iszero), Tuple{Tuple{Int64}}}})
    Base.precompile(Tuple{typeof(broadcasted),Function,ComplexF64})
    Base.precompile(Tuple{typeof(broadcasted),Function,Float64,Matrix{ComplexF64},Matrix{ComplexF64}})
    Base.precompile(Tuple{typeof(broadcasted),Function,Int64,Base.Broadcast.Broadcasted{Base.Broadcast.Style{Tuple}, Nothing, typeof(*), Tuple{Int64, Base.Broadcast.Broadcasted{Base.Broadcast.Style{Tuple}, Nothing, typeof(-), Tuple{Int64, Tuple{Int64}}}}}})
    Base.precompile(Tuple{typeof(broadcasted),Function,Int64,Base.Broadcast.Broadcasted{Base.Broadcast.Style{Tuple}, Nothing, typeof(-), Tuple{Int64, Tuple{Int64}}}})
    Base.precompile(Tuple{typeof(broadcasted),Function,Int64,Base.Broadcast.Broadcasted{Base.Broadcast.Style{Tuple}, Nothing, typeof(<<), Tuple{Int64, Base.Broadcast.Broadcasted{Base.Broadcast.Style{Tuple}, Nothing, typeof(*), Tuple{Int64, Base.Broadcast.Broadcasted{Base.Broadcast.Style{Tuple}, Nothing, typeof(-), Tuple{Int64, Tuple{Int64}}}}}}}})
    Base.precompile(Tuple{typeof(broadcasted),Function,Int64,Tuple{Int64}})
    Base.precompile(Tuple{typeof(broadcasted),Function,Matrix{ComplexF64},Matrix{ComplexF64}})
    Base.precompile(Tuple{typeof(broadcasted),Function,Matrix{ComplexF64},Matrix{Float64}})
    Base.precompile(Tuple{typeof(broadcasted),Function,Matrix{ComplexF64}})
    Base.precompile(Tuple{typeof(broadcasted),Function,Tuple{Int64}})
    Base.precompile(Tuple{typeof(collect_similar),UnitRange{Int64},Generator{UnitRange{Int64}, Fix1{Type{Symbol}, Symbol}}})
    Base.precompile(Tuple{typeof(deepcopy_internal),Vector{ComplexF64},IdDict{Any, Any}})
    Base.precompile(Tuple{typeof(dotview),Array{ComplexF64, 3},Int64,Function,Function})
    Base.precompile(Tuple{typeof(hvcat),Tuple{Int64, Int64},ComplexF64,Vararg{Number}})
    Base.precompile(Tuple{typeof(hvcat),Tuple{Int64, Int64},Float64,Vararg{Float64}})
    Base.precompile(Tuple{typeof(materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(real), Tuple{Vector{Float64}}}})
    Base.precompile(Tuple{typeof(materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(/), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(-), Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}}, Matrix{Float64}}}, Float64}}})
    Base.precompile(Tuple{typeof(materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(/), Tuple{Matrix{ComplexF64}, Matrix{Float64}}}, Float64}}})
    Base.precompile(Tuple{typeof(promote_typeof),Float64,Float64,Vararg{Any}})
    Base.precompile(Tuple{typeof(string),Symbol,Int64,Symbol,Int64})
    Base.precompile(Tuple{typeof(typed_hvcat),Type{ComplexF64},Tuple{Int64, Int64},ComplexF64,Vararg{Number}})
    Base.precompile(Tuple{typeof(vcat),Vector{Int64},Vector{Int64},Vector{Int64}})
    Base.precompile(Tuple{typeof(vcat),Vector{Int64}})
    Base.precompile(Tuple{typeof(vect),Matrix{ComplexF64}})
    Base.precompile(Tuple{typeof(|>),StepRangeLen{Float64, TwicePrecision{Float64}, TwicePrecision{Float64}, Int64},typeof(collect)})
    Base.precompile(Tuple{var"##s88#234",Any,Any,Any})
    let fbody = try __lookup_kwbody__(which(sum, (Function,Tuple{},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(sum),Function,Tuple{},))
    end
end
end
