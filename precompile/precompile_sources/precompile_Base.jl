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
    Base.precompile(Tuple{typeof(zeros),Type{ComplexF64},Union{Integer, AbstractUnitRange},Int64,Int64})   # time: 0.113732584
    Base.precompile(Tuple{typeof(hvcat),Tuple{Int64, Int64},ComplexF64,Vararg{Number}})   # time: 0.061970416
    Base.precompile(Tuple{typeof(getindex),Array{ComplexF64, 3},Int64,Int64,Colon})   # time: 0.02398879
    Base.precompile(Tuple{typeof(vcat),Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vararg{Vector}})   # time: 0.023501538
    Base.precompile(Tuple{typeof(vcat),Vararg{Vector}})   # time: 0.020803837
    Base.precompile(Tuple{typeof(repeat),Vector{Bool},UInt64})   # time: 0.01416537
    Base.precompile(Tuple{typeof(iterate),ColumnSlices{Matrix{ComplexF64}, Tuple{OneTo{Int64}}, SubArray{ComplexF64, 1, Matrix{ComplexF64}, Tuple{Slice{OneTo{Int64}}, Int64}, true}}})   # time: 0.009196581
    Base.precompile(Tuple{typeof(real),Vector})   # time: 0.007327673
    Base.precompile(Tuple{typeof(eachcol),Matrix{ComplexF64}})   # time: 0.007235681
    Base.precompile(Tuple{typeof(sum),Vector{Matrix{ComplexF64}}})   # time: 0.006800835
    Base.precompile(Tuple{typeof(setindex!),Array{ComplexF64, 3},Array{ComplexF64, 3},UnitRange{Int64},Colon,Colon})   # time: 0.005965119
    Base.precompile(Tuple{typeof(+),Vector{Float64},Vector{ComplexF64}})   # time: 0.005800171
    Base.precompile(Tuple{typeof(imag),Vector})   # time: 0.005636672
    Base.precompile(Tuple{typeof(__cat),Vector{Any},Tuple{Int64},Tuple{Bool},Vector{Vector{UInt64}},Vararg{Any}})   # time: 0.005239376
    Base.precompile(Tuple{typeof(create_expr_cache),PkgId,String,String,String,Vector{Pair{PkgId, UInt128}},IO,IO})   # time: 0.004779335
    Base.precompile(Tuple{typeof(promote_op),typeof(/),Type{Matrix{ComplexF64}},Type{Float64}})   # time: 0.004611777
    Base.precompile(Tuple{typeof(setindex_widen_up_to),Vector{Vector{Vector{UInt64}}},Vector{Any},Int64})   # time: 0.00449808
    Base.precompile(Tuple{typeof(sum),Vector{Any}})   # time: 0.004387917
    Base.precompile(Tuple{typeof(repeat),Vector{UInt64},Int64})   # time: 0.004187997
    Base.precompile(Tuple{typeof(setindex_widen_up_to),Vector{T} where T<:Complex,Complex,Int64})   # time: 0.003857708
    Base.precompile(Tuple{typeof(*),Complex{Int64},Vector{Float64}})   # time: 0.003622464
    Base.precompile(Tuple{typeof(imag),Vector{ComplexF64}})   # time: 0.003256788
    Base.precompile(Tuple{typeof(count),Fix2{typeof(!=), Int64},NTuple{6, Int64}})   # time: 0.002747001
    Base.precompile(Tuple{typeof(count),Fix2{typeof(!=), Int64},NTuple{5, Int64}})   # time: 0.002610791
    Base.precompile(Tuple{typeof(getindex),Array{ComplexF64, 3},Int64,Vector{UInt64},Vector{UInt64}})   # time: 0.00260754
    Base.precompile(Tuple{typeof(count),Fix2{typeof(!=), Int64},Tuple{Int64, Int64}})   # time: 0.002518293
    Base.precompile(Tuple{typeof(getindex),Array{ComplexF64, 3},Int64,Colon,Int64})   # time: 0.002510458
    Base.precompile(Tuple{typeof(setindex!),Vector{Vector{Any}},Vector{Vector{UInt64}},Int64})   # time: 0.002469373
    Base.precompile(Tuple{typeof(setindex_widen_up_to),Vector{Vector{Any}},Vector{Vector{UInt64}},Int64})   # time: 0.002460206
    Base.precompile(Tuple{typeof(count),Fix2{typeof(!=), Int64},NTuple{4, Int64}})   # time: 0.002331879
    Base.precompile(Tuple{typeof(count),Fix2{typeof(!=), Int64},Tuple{Int64, Int64, Int64}})   # time: 0.0022525
    Base.precompile(Tuple{typeof(==),Matrix{ComplexF64},Matrix{ComplexF64}})   # time: 0.002225163
    let fbody = try __lookup_kwbody__(which(cat, (Vector{Vector{UInt64}},Vararg{Any},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Val{1},typeof(cat),Vector{Vector{UInt64}},Vararg{Any},))
    end
end   # time: 0.002136123
    Base.precompile(Tuple{typeof(_prod),Vector{ComplexF64},Colon})   # time: 0.002061756
    Base.precompile(Tuple{typeof(setindex_widen_up_to),Vector{ComplexF64},Any,Int64})   # time: 0.002017872
    Base.precompile(Tuple{typeof(sum),Vector{ComplexF64}})   # time: 0.002014834
    Base.precompile(Tuple{typeof(append!),Vector{Any},Int64})   # time: 0.00197821
    Base.precompile(Tuple{typeof(repeat),Vector{Bool},Int64})   # time: 0.001935501
    Base.precompile(Tuple{typeof(promote_typejoin_union),Type{Union{Vector{Any}, Vector{Vector{UInt64}}}}})   # time: 0.001915657
    Base.precompile(Tuple{typeof(setindex_widen_up_to),Vector{ComplexF64},Complex,Int64})   # time: 0.001883168
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{dims::Val{1}},typeof(cat),Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vararg{Any}})   # time: 0.001856081
    Base.precompile(Tuple{typeof(any),Function,Matrix{ComplexF64}})   # time: 0.001856033
    Base.precompile(Tuple{typeof(append!),Vector{Any},Any})   # time: 0.001793499
    Base.precompile(Tuple{typeof(_prod),Vector{Int64},Colon})   # time: 0.001763455
    Base.precompile(Tuple{typeof(open),CmdRedirect,String,TTY})   # time: 0.001578221
    Base.precompile(Tuple{Colon,Int64,Real})   # time: 0.001446585
    Base.precompile(Tuple{typeof(string),ComplexF64})   # time: 0.001372203
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{digits::Int64},typeof(round),Real})   # time: 0.00126558
    Base.precompile(Tuple{typeof(prod),Vector})   # time: 0.001245206
    Base.precompile(Tuple{typeof(close),PipeEndpoint})   # time: 0.001205957
    Base.precompile(Tuple{typeof(_cat),Val{1},Vector{Vector{UInt64}},Vector{Vector{UInt64}},Vararg{Any}})   # time: 0.001175463
    Base.precompile(Tuple{typeof(_sub2ind),Tuple{OneTo{Int64}, OneTo{Int64}, OneTo{Int64}},Int64,Vararg{Integer}})   # time: 0.001044378
    Base.precompile(Tuple{typeof(vcat),Vector{Int64},Vector{Int64},Vector{Int64}})   # time: 0.001018001
end
