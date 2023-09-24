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
    Base.precompile(Tuple{Type{ITensor},Index{Int64},Vararg{Index{Int64}}})
    Base.precompile(Tuple{Type{ITensor},Type{EmptyNumber},Index{Int64},Vararg{Any}})
    Base.precompile(Tuple{Type{ITensor},Type{EmptyNumber},Tuple{Index{Int64}}})
    Base.precompile(Tuple{Type{ITensor},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:allow_alias,), Tuple{Bool}},typeof(permute),ITensor,Index{Int64},Index{Int64}})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:lefttags, :cutoff, :maxdim, :method), Tuple{TagSet, Float64, Int64, String}},typeof(svd),ITensor,Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:lefttags, :cutoff, :maxdim, :method, :nsite, :nsweeps), Tuple{TagSet, Float64, Int64, String, Int64, Int64}},typeof(svd),ITensor,Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:ortho, :which_decomp), Tuple{String, String}},typeof(factorize),ITensor,Index{Int64}})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:tags,), Tuple{TagSet}},typeof(factorize),ITensor,Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:tags,), Tuple{TagSet}},typeof(qr),ITensor,Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Base.ReshapedArray{ComplexF64, 1, Adjoint{ComplexF64, Matrix{ComplexF64}}, Tuple{Base.MultiplicativeInverses.SignedMultiplicativeInverse{Int64}}}}},Tensor{Number, 2, Combiner, Tuple{Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Base.ReshapedArray{ComplexF64, 1, Adjoint{ComplexF64, Matrix{ComplexF64}}, Tuple{Base.MultiplicativeInverses.SignedMultiplicativeInverse{Int64}}}}},Tensor{Number, 3, Combiner, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{Float64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{Float64, Vector{Float64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tensor{Number, 2, Combiner, Tuple{Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tensor{Number, 3, Combiner, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tensor{Number, 2, Combiner, Tuple{Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tensor{Number, 3, Combiner, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{Int64, 1, Tuple{Index{Int64}}, Dense{Int64, Vector{Int64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tensor{Number, 3, Combiner, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 5, NTuple{5, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 5, NTuple{5, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 5, NTuple{5, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{Int64, 1, Tuple{Index{Int64}}, Dense{Int64, Vector{Int64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 5, NTuple{5, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tensor{Number, 3, Combiner, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 6, NTuple{6, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 6, NTuple{6, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{Float64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{Float64, Vector{Float64}}},DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DiagTensor{Float64, 2, Tuple{Index{Int64}, Index{Int64}}, Diag{Float64, Vector{Float64}}},DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DiagTensor{Float64, 2, Tuple{Index{Int64}, Index{Int64}}, Diag{Float64, Vector{Float64}}},DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_getindex),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_intersect),NTuple{4, Index{Int64}},NTuple{4, Index{Int64}}})
    Base.precompile(Tuple{typeof(_intersect),NTuple{4, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_intersect),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},NTuple{4, Index{Int64}}})
    Base.precompile(Tuple{typeof(_intersect),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_intersect),Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_intersect),Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}}})
    Base.precompile(Tuple{typeof(_intersect),Tuple{Index{Int64}, Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(_permute),AllowAlias,DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_permute),NeverAlias,DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff!),Vector{Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff!),Vector{Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff!),Vector{Index{Int64}},Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Vararg{Any}})
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}},Vararg{Any}})
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}},Vararg{Tuple{Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, UInt8}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Int64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, UInt8}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{Float64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{Float64, Vector{Float64}}},Float64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{Int64, 1, Tuple{Index{Int64}}, Dense{Int64, Vector{Int64}}},Int64,Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),EmptyTensor{EmptyNumber, 1, Tuple{Index{Int64}}, EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}}},Int64,Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),EmptyTensor{EmptyNumber, 2, Tuple{Index{Int64}, Index{Int64}}, EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),EmptyTensor{EmptyNumber, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),EmptyTensor{EmptyNumber, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, UInt8}})
    Base.precompile(Tuple{typeof(_setindex!!),EmptyTensor{EmptyNumber, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}}},Float64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),EmptyTensor{EmptyNumber, 4, NTuple{4, Index{Int64}}, EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),EmptyTensor{EmptyNumber, 4, NTuple{4, Index{Int64}}, EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, UInt8}})
    Base.precompile(Tuple{typeof(_union),Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(combiner),Index{Int64},Index{Int64}})
    Base.precompile(Tuple{typeof(combiner),Index{Int64}})
    Base.precompile(Tuple{typeof(combiner),Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(dag),AllowAlias,Tensor{Number, 2, Combiner, Tuple{Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(dag),AllowAlias,Tensor{Number, 3, Combiner, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(eachindval),Index{Int64}})
    Base.precompile(Tuple{typeof(filter_inds_set_function),Function,Function,ITensor,ITensor,ITensor})
    Base.precompile(Tuple{typeof(filter_inds_set_function),Function,Function,ITensor,Tuple{Index{Int64}}})
    Base.precompile(Tuple{typeof(filter_inds_set_function),Function,ITensor,ITensor,ITensor,Vararg{ITensor}})
    Base.precompile(Tuple{typeof(filter_inds_set_function),Function,ITensor,ITensor})
    Base.precompile(Tuple{typeof(indices),Index{Int64},Index{Int64},Vararg{Index{Int64}}})
    Base.precompile(Tuple{typeof(noprime),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(permute),ITensor,Index{Int64},Index{Int64}})
    Base.precompile(Tuple{typeof(prime),DenseTensor{Int64, 1, Tuple{Index{Int64}}, Dense{Int64, Vector{Int64}}}})
    Base.precompile(Tuple{typeof(replaceind!),ITensor,Index{Int64},Index{Int64}})
    Base.precompile(Tuple{typeof(replaceind),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Index{Int64},Index{Int64}})
    Base.precompile(Tuple{typeof(replaceind),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Index{Int64},Index{Int64}})
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),ITensor,Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),NTuple{4, Index{Int64}},NTuple{4, Index{Int64}},NTuple{4, Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},NTuple{4, Index{Int64}},NTuple{4, Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),Tuple{Index{Int64}},Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceprime),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Pair{Int64, Int64}})
    Base.precompile(Tuple{typeof(replaceprime),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Pair{Int64, Int64}})
    Base.precompile(Tuple{typeof(setindex!),ITensor,Int64,Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(settags),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},TagSet,Index{Int64}})
    Base.precompile(Tuple{typeof(settags),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},TagSet,Index{Int64}})
    Base.precompile(Tuple{typeof(settags),DiagTensor{Float64, 2, Tuple{Index{Int64}, Index{Int64}}, Diag{Float64, Vector{Float64}}},TagSet,Index{Int64}})
    Base.precompile(Tuple{typeof(swapinds),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(swapinds),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(swapinds),DenseTensor{Int64, 1, Tuple{Index{Int64}}, Dense{Int64, Vector{Int64}}},Vector{Index{Int64}},Vector{Index{Int64}}})
    isdefined(ITensors, Symbol("#494#509")) && Base.precompile(Tuple{getfield(ITensors, Symbol("#494#509")),Vector{Index{Int64}}})
    let fbody = try __lookup_kwbody__(which(commonind, (ITensor,Vararg{ITensor},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(commonind),ITensor,Vararg{ITensor},))
    end
end
    let fbody = try __lookup_kwbody__(which(filter_inds_set_function, (Function,ITensor,Vararg{Any},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(filter_inds_set_function),Function,ITensor,Vararg{Any},))
    end
end
    let fbody = try __lookup_kwbody__(which(orthogonalize!, (MPS,Int64,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(orthogonalize!),MPS,Int64,))
    end
end
end
