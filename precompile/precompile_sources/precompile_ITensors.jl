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
    Base.precompile(Tuple{Type{ITensor},AllowAlias,Combiner,Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{Type{ITensor},AllowAlias,DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{Type{ITensor},AllowAlias,EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}},Tuple{Index{Int64}}})
    Base.precompile(Tuple{Type{ITensor},Index{Int64},Vararg{Index{Int64}}})
    Base.precompile(Tuple{Type{ITensor},Type{EmptyNumber},Index{Int64},Vararg{Any}})
    Base.precompile(Tuple{Type{ITensor},Vector{Index{Int64}}})
    Base.precompile(Tuple{Type{TruncEigen},ITensor,ITensor,ITensor,Spectrum{Vector{Float64}, Float64},Index{Int64},Index{Int64}})
    Base.precompile(Tuple{Type{TruncSVD},ITensor,ITensor,ITensor,Spectrum{Vector{Float64}, Float64},Index{Int64},Index{Int64}})
    Base.precompile(Tuple{Type{Vector{IndexT} where IndexT<:Index},Index{Int64},Vararg{Index{Int64}}})
    Base.precompile(Tuple{typeof(*),Bool,ITensor})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{allow_alias::Bool},typeof(permute),ITensor,Index{Int64},Index{Int64}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{cutoff::Float64, maxdim::Int64, alg::String},typeof(product),MPO,MPS})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{dir::Arrow, tags::String},typeof(combiner),Index{Int64},Index{Int64}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{dir::Arrow, tags::String},typeof(combiner),Index{Int64}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{dir::Nothing},typeof(combiner),Index{Int64},Index{Int64}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{dir::Nothing},typeof(combiner),Index{Int64}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{ishermitian::Bool, mindim::Int64, maxdim::Int64, cutoff::Float64, tags::TagSet, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing},typeof(eigen),ITensor,Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{ishermitian::Bool, tags::TagSet, cutoff::Float64, maxdim::Int64, mindim::Int64},typeof(eigen),ITensor,Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{lefttags::TagSet, cutoff::Float64, maxdim::Int64},typeof(svd),ITensor,Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{maxdim::Int64, mindim::Int64, cutoff::Float64, eigen_perturbation::Nothing, ortho::String, normalize::Bool, which_decomp::Nothing, svd_alg::Nothing},typeof(replacebond!),MPS,Int64,ITensor})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{mindim::Int64, maxdim::Int64, cutoff::Float64, ortho::String, which_decomp::Nothing, eigen_perturbation::Nothing, svd_alg::Nothing, tags::TagSet, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing, min_blockdim::Nothing},typeof(factorize),ITensor,Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{mindim::Int64, maxdim::Int64, cutoff::Float64, ortho::String, which_decomp::Nothing, eigen_perturbation::Nothing, svd_alg::Nothing, tags::TagSet, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing, min_blockdim::Nothing},typeof(factorize),ITensor,Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{mindim::Int64, maxdim::Int64, cutoff::Float64, tags::TagSet, ortho::String, eigen_perturbation::Nothing, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing},typeof(factorize_eigen),ITensor,Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{mindim::Int64, maxdim::Int64, cutoff::Float64, tags::TagSet, ortho::String, eigen_perturbation::Nothing, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing},typeof(factorize_eigen),ITensor,Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{mindim::Nothing, maxdim::Int64, cutoff::Float64, tags::TagSet, ortho::String, alg::Nothing, dir::Nothing, singular_values!::Nothing, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing, min_blockdim::Nothing},typeof(factorize_svd),ITensor,Index{Int64}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{mindim::Nothing, maxdim::Int64, cutoff::Nothing, tags::TagSet, ortho::String, alg::Nothing, dir::Nothing, singular_values!::Nothing, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing, min_blockdim::Nothing},typeof(factorize_svd),ITensor,Index{Int64}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{ortho::String, which_decomp::String},typeof(factorize),ITensor,Index{Int64}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{tags::TagSet, maxdim::Nothing},typeof(factorize),ITensor,Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{tags::TagSet, positive::Bool},typeof(qr),ITensor,Vector{Index{Int64}},Vector{Index{Int64}}})
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
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{Float64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{Float64, Vector{Float64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tensor{Number, 2, Combiner, Tuple{Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tensor{Number, 3, Combiner, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(_contract),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},DenseTensor{Float64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{Float64, Vector{Float64}}}})
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
    Base.precompile(Tuple{typeof(_intersect),NTuple{4, Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(_intersect),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},NTuple{4, Index{Int64}}})
    Base.precompile(Tuple{typeof(_intersect),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_intersect),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_intersect),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}}})
    Base.precompile(Tuple{typeof(_intersect),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(_intersect),Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_intersect),Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_intersect),Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}}})
    Base.precompile(Tuple{typeof(_intersect),Tuple{Index{Int64}, Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(_permute),AllowAlias,DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_permute),NeverAlias,DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff!),Vector{Index{Int64}},NTuple{4, Index{Int64}},NTuple{4, Index{Int64}},Vararg{Any}})
    Base.precompile(Tuple{typeof(_setdiff!),Vector{Index{Int64}},NTuple{4, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Vararg{Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(_setdiff!),Vector{Index{Int64}},NTuple{4, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff!),Vector{Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}},NTuple{4, Index{Int64}},Vararg{Any}})
    Base.precompile(Tuple{typeof(_setdiff!),Vector{Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Vararg{Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(_setdiff!),Vector{Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff!),Vector{Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff!),Vector{Index{Int64}},Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff!),Vector{Index{Int64}},Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff),NTuple{4, Index{Int64}},NTuple{4, Index{Int64}},Vararg{Any}})
    Base.precompile(Tuple{typeof(_setdiff),NTuple{4, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Vararg{Any}})
    Base.precompile(Tuple{typeof(_setdiff),NTuple{4, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Vararg{Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(_setdiff),NTuple{4, Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},NTuple{4, Index{Int64}},Vararg{Any}})
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Vararg{Any}})
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}},Vararg{Any}})
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}},Vararg{Tuple{Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(_setdiff),Tuple{Index{Int64}, Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, UInt64}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Int64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, UInt64}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{Float64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{Float64, Vector{Float64}}},Float64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),DenseTensor{Int64, 1, Tuple{Index{Int64}}, Dense{Int64, Vector{Int64}}},Int64,Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),EmptyTensor{EmptyNumber, 1, Tuple{Index{Int64}}, EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}}},Int64,Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),EmptyTensor{EmptyNumber, 2, Tuple{Index{Int64}, Index{Int64}}, EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),EmptyTensor{EmptyNumber, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),EmptyTensor{EmptyNumber, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, UInt64}})
    Base.precompile(Tuple{typeof(_setindex!!),EmptyTensor{EmptyNumber, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}}},Float64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),EmptyTensor{EmptyNumber, 4, NTuple{4, Index{Int64}}, EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(_setindex!!),EmptyTensor{EmptyNumber, 4, NTuple{4, Index{Int64}}, EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}}},ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, UInt64}})
    Base.precompile(Tuple{typeof(_union),Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(add_trivial_index),ITensor,Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(adjoint),ITensor})
    Base.precompile(Tuple{typeof(combiner),Index{Int64},Index{Int64}})
    Base.precompile(Tuple{typeof(combiner),Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(commonind),ITensor,Vararg{ITensor}})
    Base.precompile(Tuple{typeof(commoninds),ITensor,Vararg{Any}})
    Base.precompile(Tuple{typeof(dag),AllowAlias,DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(dag),AllowAlias,DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(dag),AllowAlias,DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(dag),AllowAlias,DenseTensor{Float64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{Float64, Vector{Float64}}}})
    Base.precompile(Tuple{typeof(dag),AllowAlias,Tensor{Number, 2, Combiner, Tuple{Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(dag),AllowAlias,Tensor{Number, 3, Combiner, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(dag),Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(eachindval),Index{Int64}})
    Base.precompile(Tuple{typeof(filter_inds_set_function),Function,Function,ITensor,ITensor,ITensor,ITensor})
    Base.precompile(Tuple{typeof(filter_inds_set_function),Function,Function,ITensor,ITensor,ITensor})
    Base.precompile(Tuple{typeof(filter_inds_set_function),Function,Function,ITensor,ITensor})
    Base.precompile(Tuple{typeof(filter_inds_set_function),Function,Function,ITensor,Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(filter_inds_set_function),Function,Function,ITensor,Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(filter_inds_set_function),Function,Function,ITensor,Tuple{Index{Int64}}})
    Base.precompile(Tuple{typeof(filter_inds_set_function),Function,Function,ITensor,Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(filter_inds_set_function),Function,ITensor,ITensor,ITensor,Vararg{ITensor}})
    Base.precompile(Tuple{typeof(filter_inds_set_function),Function,ITensor,ITensor})
    Base.precompile(Tuple{typeof(getfirst),Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(hasind),ITensor,Index{Int64}})
    Base.precompile(Tuple{typeof(indices),Index{Int64},Index{Int64},Vararg{Index{Int64}}})
    Base.precompile(Tuple{typeof(itensor),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Base.ReshapedArray{ComplexF64, 1, Adjoint{ComplexF64, Matrix{ComplexF64}}, Tuple{Base.MultiplicativeInverses.SignedMultiplicativeInverse{Int64}}}}}})
    Base.precompile(Tuple{typeof(itensor),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(itensor),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(itensor),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(itensor),DenseTensor{ComplexF64, 5, NTuple{5, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(itensor),DenseTensor{ComplexF64, 6, NTuple{6, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(itensor),DenseTensor{Float64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{Float64, Vector{Float64}}}})
    Base.precompile(Tuple{typeof(itensor),DenseTensor{Int64, 1, Tuple{Index{Int64}}, Dense{Int64, Vector{Int64}}}})
    Base.precompile(Tuple{typeof(itensor),DiagTensor{Float64, 2, Tuple{Index{Int64}, Index{Int64}}, Diag{Float64, Vector{Float64}}}})
    Base.precompile(Tuple{typeof(itensor),Tensor{Number, 2, Combiner, Tuple{Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(itensor),Tensor{Number, 3, Combiner, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(noprime),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}}})
    Base.precompile(Tuple{typeof(noprime),ITensor})
    Base.precompile(Tuple{typeof(permute),ITensor,Index{Int64},Index{Int64}})
    Base.precompile(Tuple{typeof(permute),Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(prime),DenseTensor{Int64, 1, Tuple{Index{Int64}}, Dense{Int64, Vector{Int64}}}})
    Base.precompile(Tuple{typeof(remove_trivial_index),ITensor,ITensor,Nothing,Nothing})
    Base.precompile(Tuple{typeof(replaceind!),ITensor,Index{Int64},Index{Int64}})
    Base.precompile(Tuple{typeof(replaceind),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Index{Int64},Index{Int64}})
    Base.precompile(Tuple{typeof(replaceind),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Index{Int64},Index{Int64}})
    Base.precompile(Tuple{typeof(replaceind),ITensor,Index{Int64},Index{Int64}})
    Base.precompile(Tuple{typeof(replaceinds!),ITensor,Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Pair{Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Pair{Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Pair{Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{Float64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{Float64, Vector{Float64}}},Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),DenseTensor{Float64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{Float64, Vector{Float64}}},Vector{Pair{Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(replaceinds),DiagTensor{Float64, 2, Tuple{Index{Int64}, Index{Int64}}, Diag{Float64, Vector{Float64}}},Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),ITensor,Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),ITensor,Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),ITensor,Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),ITensor,Vector{Pair{Index{Int64}, Index{Int64}}}})
    Base.precompile(Tuple{typeof(replaceinds),NTuple{4, Index{Int64}},NTuple{4, Index{Int64}},NTuple{4, Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),NTuple{4, Index{Int64}},Pair{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},NTuple{4, Index{Int64}},NTuple{4, Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),Tuple{Index{Int64}, Index{Int64}, Index{Int64}},Pair{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),Tuple{Index{Int64}, Index{Int64}},Pair{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceinds),Tuple{Index{Int64}},Tuple{Index{Int64}, Index{Int64}},Tuple{Index{Int64}, Index{Int64}}})
    Base.precompile(Tuple{typeof(replaceprime),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Pair{Int64, Int64}})
    Base.precompile(Tuple{typeof(replaceprime),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Pair{Int64, Int64}})
    Base.precompile(Tuple{typeof(setindex!),ITensor,ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(setindex!),ITensor,ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, UInt64}})
    Base.precompile(Tuple{typeof(setindex!),ITensor,ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(setindex!),ITensor,ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, UInt64}})
    Base.precompile(Tuple{typeof(setindex!),ITensor,ComplexF64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(setindex!),ITensor,Int64,Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64},Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(setindex!),ITensor,Int64,Pair{Index{Int64}, Int64}})
    Base.precompile(Tuple{typeof(settags),DenseTensor{ComplexF64, 2, Tuple{Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},TagSet,Index{Int64}})
    Base.precompile(Tuple{typeof(settags),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},TagSet,Index{Int64}})
    Base.precompile(Tuple{typeof(settags),DiagTensor{Float64, 2, Tuple{Index{Int64}, Index{Int64}}, Diag{Float64, Vector{Float64}}},TagSet,Index{Int64}})
    Base.precompile(Tuple{typeof(settags),ITensor,TagSet,Index{Int64}})
    Base.precompile(Tuple{typeof(swapinds!),ITensor,Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(swapinds),DenseTensor{ComplexF64, 3, Tuple{Index{Int64}, Index{Int64}, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(swapinds),DenseTensor{ComplexF64, 4, NTuple{4, Index{Int64}}, Dense{ComplexF64, Vector{ComplexF64}}},Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(swapinds),DenseTensor{Int64, 1, Tuple{Index{Int64}}, Dense{Int64, Vector{Int64}}},Vector{Index{Int64}},Vector{Index{Int64}}})
    Base.precompile(Tuple{typeof(uniqueinds),ITensor,Vararg{Any}})
    isdefined(ITensors, Symbol("#42#43")) && Base.precompile(Tuple{getfield(ITensors, Symbol("#42#43")),Int64})
    isdefined(ITensors, Symbol("#496#511")) && Base.precompile(Tuple{getfield(ITensors, Symbol("#496#511")),Vector{Index{Int64}}})
    isdefined(ITensors, Symbol("#498#512")) && Base.precompile(Tuple{getfield(ITensors, Symbol("#498#512")),Vector{Index{Int64}}})
    isdefined(ITensors, Symbol("#506#516")) && Base.precompile(Tuple{getfield(ITensors, Symbol("#506#516")),Vector{Index{Int64}}})
    isdefined(ITensors, Symbol("#508#517")) && Base.precompile(Tuple{getfield(ITensors, Symbol("#508#517")),Vector{Index{Int64}}})
    let fbody = try __lookup_kwbody__(which(commonind, (ITensor,Vararg{Any},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, @NamedTuple{}},typeof(commonind),ITensor,Vararg{Any},))
    end
end
    let fbody = try __lookup_kwbody__(which(commonind, (ITensor,Vararg{ITensor},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, @NamedTuple{}},typeof(commonind),ITensor,Vararg{ITensor},))
    end
end
    let fbody = try __lookup_kwbody__(which(filter_inds_set_function, (Function,ITensor,Vararg{Any},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, @NamedTuple{}},typeof(filter_inds_set_function),Function,ITensor,Vararg{Any},))
    end
end
    let fbody = try __lookup_kwbody__(which(swapinds, (ITensor,Vector{Index{Int64}},Vararg{Vector{Index{Int64}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, @NamedTuple{}},typeof(swapinds),ITensor,Vector{Index{Int64}},Vararg{Vector{Index{Int64}}},))
    end
end
    let fbody = try __lookup_kwbody__(which(uniqueind, (ITensor,Vararg{Any},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, @NamedTuple{}},typeof(uniqueind),ITensor,Vararg{Any},))
    end
end
    let fbody = try __lookup_kwbody__(which(uniqueind, (ITensor,Vararg{ITensor},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, @NamedTuple{}},typeof(uniqueind),ITensor,Vararg{ITensor},))
    end
end
end
