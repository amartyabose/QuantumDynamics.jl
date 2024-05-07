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
    Base.precompile(Tuple{Type{CartesianIndices},Tuple{Int64}})
    Base.precompile(Tuple{Type{Fix1},Type{MappingRF},Type})
    Base.precompile(Tuple{Type{Generator},Fix1{Type{Symbol}, Symbol},UnitRange{Int64}})
    Base.precompile(Tuple{typeof(!=),NTuple{4, Int64},NTuple{4, Int64}})
    Base.precompile(Tuple{typeof(*),ComplexF64,ComplexF64,ComplexF64,ComplexF64})
    Base.precompile(Tuple{typeof(*),Complex{Int64},Vector{Float64}})
    Base.precompile(Tuple{typeof(*),Float64,Float16})
    Base.precompile(Tuple{typeof(+),Vector{Float64},Vector{ComplexF64}})
    Base.precompile(Tuple{typeof(==),Matrix{ComplexF64},Matrix{ComplexF64}})
    Base.precompile(Tuple{typeof(Base.Cartesian._nloops),Int64,Symbol,Symbol,Expr,Vararg{Expr}})
    Base.precompile(Tuple{typeof(Base.MainInclude.include),String})
    Base.precompile(Tuple{typeof(Base.SimdLoop.compile),Expr,Symbol})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{dims::Val{1}},typeof(cat),Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Any},Vector{Any},Vector{Any},Vector{Any}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{dims::Val{1}},typeof(cat),Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Any},Vector{Any},Vector{Any}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{dims::Val{1}},typeof(cat),Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Any},Vector{Any}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{dims::Val{1}},typeof(cat),Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Any}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{init::Int64},typeof(mapreduce),Type,Function,Tuple{Int64}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{rev::Bool, by::typeof(abs)},typeof(sortperm),Vector{Float64}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{stderr::TTY, stdout::TTY},typeof(pipeline),Cmd})
    Base.precompile(Tuple{typeof(__cat),Vector{Any},Tuple{Int64},Tuple{Bool},Vector{Vector{UInt64}},Vararg{Any}})
    Base.precompile(Tuple{typeof(__cat_offset!),Vector{Any},Tuple{Int64},Tuple{Bool},Tuple{Int64},Vector{Any},Vector{Any},Vararg{Vector{Any}}})
    Base.precompile(Tuple{typeof(__cat_offset!),Vector{Any},Tuple{Int64},Tuple{Bool},Tuple{Int64},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vararg{Any}})
    Base.precompile(Tuple{typeof(_array_for),Type{Tuple{Int64}},HasLength,Int64})
    Base.precompile(Tuple{typeof(_cat_size_shape),Tuple{Bool},Tuple{Int64},Vector{Any},Vector{Any},Vararg{Vector{Any}}})
    Base.precompile(Tuple{typeof(_cat_size_shape),Tuple{Bool},Tuple{Int64},Vector{Vector{UInt64}},Vector{Any},Vararg{Vector{Any}}})
    Base.precompile(Tuple{typeof(_cat_size_shape),Tuple{Bool},Tuple{Int64},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vararg{Any}})
    Base.precompile(Tuple{typeof(_cat_t),Val{1},Type{Any},Vector{Vector{UInt64}},Vararg{Any}})
    Base.precompile(Tuple{typeof(_nt_names),Type{@NamedTuple{cutoff::Float64, maxdim::Int64}}})
    Base.precompile(Tuple{typeof(_nt_names),Type{@NamedTuple{mindim::Int64, maxdim::Int64, cutoff::Float64, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing}}})
    Base.precompile(Tuple{typeof(_nt_names),Type{@NamedTuple{positive::Bool}}})
    Base.precompile(Tuple{typeof(_nt_names),Type{@NamedTuple{reltol::Float64, abstol::Float64, saveat::Float64}}})
    Base.precompile(Tuple{typeof(_nt_names),Type{@NamedTuple{wait::Bool}}})
    Base.precompile(Tuple{typeof(_prod),Vector{Int64},Colon})
    Base.precompile(Tuple{typeof(_tryrequire_from_serialized),PkgId,String,String})
    Base.precompile(Tuple{typeof(all),Function,SimpleVector})
    Base.precompile(Tuple{typeof(any),Function,Matrix{ComplexF64}})
    Base.precompile(Tuple{typeof(append!),Vector{Any},Vector{Symbol}})
    Base.precompile(Tuple{typeof(argtail),Type})
    Base.precompile(Tuple{typeof(broadcasted),typeof(*),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Nothing, typeof(*), Tuple{Complex{Int64}, Float64}},Vector{Float64}})
    Base.precompile(Tuple{typeof(broadcasted),typeof(*),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(/), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(-), Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}}, Matrix{Float64}}},Float64})
    Base.precompile(Tuple{typeof(broadcasted),typeof(*),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(/), Tuple{Matrix{ComplexF64}, Matrix{Float64}}},Float64})
    Base.precompile(Tuple{typeof(broadcasted),typeof(*),Complex{Int64},Float64})
    Base.precompile(Tuple{typeof(broadcasted),typeof(*),Vector{ComplexF64},Vector{ComplexF64}})
    Base.precompile(Tuple{typeof(broadcasted),typeof(*),Vector{ComplexF64},Vector{Float64}})
    Base.precompile(Tuple{typeof(broadcasted),typeof(*),Vector{Float64},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Float64, Vector{Float64}}}, Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Nothing, typeof(*), Tuple{Complex{Int64}, Float64}}, Vector{Float64}}}}}})
    Base.precompile(Tuple{typeof(broadcasted),typeof(+),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Float64, Vector{Float64}}},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Nothing, typeof(*), Tuple{Complex{Int64}, Float64}}, Vector{Float64}}}})
    Base.precompile(Tuple{typeof(broadcasted),typeof(+),Vector{ComplexF64},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Vector{Float64}, Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Float64, Vector{Float64}}}, Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Nothing, typeof(*), Tuple{Complex{Int64}, Float64}}, Vector{Float64}}}}}}}})
    Base.precompile(Tuple{typeof(broadcasted),typeof(/),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(-), Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}},Matrix{Float64}})
    Base.precompile(Tuple{typeof(broadcasted),typeof(exp),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Float64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(broadcasted),typeof(exp),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Vector{Float64}, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(broadcasted),typeof(exp),ComplexF64})
    Base.precompile(Tuple{typeof(broadcasted),typeof(muladd),Float64,Matrix{ComplexF64},Matrix{ComplexF64}})
    Base.precompile(Tuple{typeof(close),PipeEndpoint})
    Base.precompile(Tuple{typeof(collect_similar),UnitRange{Int64},Generator{UnitRange{Int64}, Fix1{Type{Symbol}, Symbol}}})
    Base.precompile(Tuple{typeof(create_expr_cache),PkgId,String,String,String,Vector{Pair{PkgId, UInt128}},IO,IO})
    Base.precompile(Tuple{typeof(deepcopy_internal),Vector{ComplexF64},IdDict{Any, Any}})
    Base.precompile(Tuple{typeof(dotview),Array{ComplexF64, 3},Int64,Function,Function})
    Base.precompile(Tuple{typeof(eltype),Type{Union{}}})
    Base.precompile(Tuple{typeof(enumerate),CartesianIndices{1, Tuple{OneTo{Int64}}}})
    Base.precompile(Tuple{typeof(enumerate),Vector{Tuple{Int64}}})
    Base.precompile(Tuple{typeof(foreach),typeof(invokelatest),Vector{Function}})
    Base.precompile(Tuple{typeof(getindex),Array{ComplexF64, 3},Int64,Function,Function})
    Base.precompile(Tuple{typeof(getindex),Array{ComplexF64, 3},Int64,Vector{UInt64},Vector{UInt64}})
    Base.precompile(Tuple{typeof(getindex),Dict{PkgId, Vector{Function}},PkgId})
    Base.precompile(Tuple{typeof(getindex),Matrix{ComplexF64},Function,UnitRange{Int64}})
    Base.precompile(Tuple{typeof(getindex),Matrix{ComplexF64},Function,Vector{Int64}})
    Base.precompile(Tuple{typeof(getindex),Tuple,UnitRange{Int64}})
    Base.precompile(Tuple{typeof(getindex),Vector{Vector{UInt64}},UInt64})
    Base.precompile(Tuple{typeof(hvcat),Tuple{Int64, Int64},ComplexF64,Vararg{Number}})
    Base.precompile(Tuple{typeof(hvcat),Tuple{Int64, Int64},Float64,Vararg{Float64}})
    Base.precompile(Tuple{typeof(hvcat_fill!),Matrix{ComplexF64},Tuple{ComplexF64, Float64, Float64, Float64}})
    Base.precompile(Tuple{typeof(imag),Vector{ComplexF64}})
    Base.precompile(Tuple{typeof(indexed_iterate),Tuple{StepRangeLen{Float64, TwicePrecision{Float64}, TwicePrecision{Float64}, Int64}, Array{ComplexF64, 3}},Int64})
    Base.precompile(Tuple{typeof(indexed_iterate),Tuple{Vector{Float64}, Array{ComplexF64, 3}},Int64})
    Base.precompile(Tuple{typeof(iterate),Enumerate{CartesianIndices{1, Tuple{OneTo{Int64}}}},Tuple{Int64, CartesianIndex{1}}})
    Base.precompile(Tuple{typeof(iterate),Enumerate{CartesianIndices{1, Tuple{OneTo{Int64}}}},Tuple{Int64}})
    Base.precompile(Tuple{typeof(map),Type,Nothing})
    Base.precompile(Tuple{typeof(materialize!),SubArray{ComplexF64, 1, Array{ComplexF64, 3}, Tuple{Int64, Int64, Slice{OneTo{Int64}}}, true},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Vector{ComplexF64}, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(materialize!),SubArray{ComplexF64, 1, Array{ComplexF64, 3}, Tuple{Int64, Slice{OneTo{Int64}}, Int64}, true},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Vector{ComplexF64}, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(materialize!),SubArray{ComplexF64, 2, Array{ComplexF64, 3}, Tuple{Int64, Vector{UInt64}, Vector{UInt64}}, false},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(identity), Tuple{Matrix{ComplexF64}}}})
    Base.precompile(Tuple{typeof(materialize!),Vector{ComplexF64},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(+), Tuple{Vector{ComplexF64}, Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Vector{Float64}, Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Float64, Vector{Float64}}}, Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Nothing, typeof(*), Tuple{Complex{Int64}, Float64}}, Vector{Float64}}}}}}}}}})
    Base.precompile(Tuple{typeof(materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Nothing, typeof(exp), Tuple{ComplexF64}}})
    Base.precompile(Tuple{typeof(materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Vector{ComplexF64}, Vector{Float64}}}})
    Base.precompile(Tuple{typeof(materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(exp), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Float64, Vector{ComplexF64}}}}}})
    Base.precompile(Tuple{typeof(materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(exp), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(*), Tuple{Vector{Float64}, Vector{ComplexF64}}}}}})
    Base.precompile(Tuple{typeof(materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(real), Tuple{Vector{Float64}}}})
    Base.precompile(Tuple{typeof(materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(/), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(-), Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}}, Matrix{Float64}}}, Float64}}})
    Base.precompile(Tuple{typeof(materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(/), Tuple{Matrix{ComplexF64}, Matrix{Float64}}}, Float64}}})
    Base.precompile(Tuple{typeof(maybeview),Array{ComplexF64, 3},Int64,Vector{UInt64},Vector{UInt64}})
    Base.precompile(Tuple{typeof(open),CmdRedirect,String,TTY})
    Base.precompile(Tuple{typeof(promote_eltypeof),Vector{Vector{UInt64}},Vector{Any},Vararg{Vector{Any}}})
    Base.precompile(Tuple{typeof(promote_typeof),Float64,Float64,Vararg{Any}})
    Base.precompile(Tuple{typeof(setindex!),Array{ComplexF64, 3},ComplexF64,Int64,UInt64,UInt64})
    Base.precompile(Tuple{typeof(similar),CartesianIndices{1, Tuple{OneTo{Int64}}},Type{Expr}})
    Base.precompile(Tuple{typeof(spawn_opts_inherit),DevNull,TTY,TTY})
    Base.precompile(Tuple{typeof(string),Symbol,Int64,Symbol,Int64})
    Base.precompile(Tuple{typeof(sum),Function,Tuple{}})
    Base.precompile(Tuple{typeof(sum),Vector{ComplexF64}})
    Base.precompile(Tuple{typeof(sym_in),Symbol,NTuple{14, Symbol}})
    Base.precompile(Tuple{typeof(sym_in),Symbol,NTuple{16, Symbol}})
    Base.precompile(Tuple{typeof(sym_in),Symbol,NTuple{51, Symbol}})
    Base.precompile(Tuple{typeof(vcat),Vector{Int64},Vector{Int64},Vector{Int64}})
    Base.precompile(Tuple{typeof(vcat),Vector{Int64}})
    Base.precompile(Tuple{typeof(vcat),Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vararg{Vector{Vector{UInt64}}}})
    Base.precompile(Tuple{typeof(vcat),Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vararg{Vector}})
    Base.precompile(Tuple{typeof(vcat),Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Vector{UInt64}}})
    Base.precompile(Tuple{typeof(vcat),Vector{Vector{UInt64}},Vector{Vector{UInt64}}})
    Base.precompile(Tuple{typeof(vcat),Vector{Vector{UInt64}}})
    Base.precompile(Tuple{typeof(vect),Matrix{ComplexF64}})
    Base.precompile(Tuple{typeof(write),PipeEndpoint,String})
    let fbody = try __lookup_kwbody__(which(sum, (Function,Tuple{},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Pairs{Symbol, Union{}, Tuple{}, @NamedTuple{}},typeof(sum),Function,Tuple{},))
    end
end
end
