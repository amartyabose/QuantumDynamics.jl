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
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{mindim::Int64, maxdim::Int64, cutoff::Float64, use_absolute_cutoff::Nothing, use_relative_cutoff::Nothing},typeof(truncate!!),Vector{Float64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{2}},Tuple{Int64, Int64, Int64},Tuple{Int64, Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{2}},Tuple{Int64, Int64},Tuple{Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{3}},NTuple{4, Int64},Tuple{Int64, Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{3}},NTuple{4, Int64},Tuple{Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{3}},Tuple{Int64, Int64, Int64},NTuple{4, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{3}},Tuple{Int64, Int64, Int64},Tuple{Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{3}},Tuple{Int64, Int64},Tuple{Int64, Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{4}},NTuple{4, Int64},NTuple{4, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{4}},NTuple{4, Int64},Tuple{Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{4}},NTuple{5, Int64},Tuple{Int64, Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{4}},NTuple{5, Int64},Tuple{Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{4}},Tuple{Int64, Int64, Int64},NTuple{5, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{4}},Tuple{Int64, Int64, Int64},Tuple{Int64, Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{5}},NTuple{4, Int64},Tuple{Int64, Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{5}},NTuple{5, Int64},NTuple{4, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{5}},NTuple{6, Int64},Tuple{Int64, Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{5}},Tuple{Int64, Int64, Int64},NTuple{4, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{6}},NTuple{5, Int64},Tuple{Int64, Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{6}},NTuple{6, Int64},NTuple{4, Int64}})
    Base.precompile(Tuple{typeof(getindex),DenseTensor{Bool, 0, Tuple{}, Dense{Bool, Vector{Bool}}}})
    Base.precompile(Tuple{typeof(intersect_positions),Tuple{Int64, Int64, Int64},NTuple{4, Int64}})
    Base.precompile(Tuple{typeof(intersect_positions),Tuple{Int64, Int64, Int64},Tuple{Int64, Int64}})
    Base.precompile(Tuple{typeof(mul!!),Matrix{ComplexF64},Matrix{ComplexF64},Transpose{ComplexF64, Matrix{ComplexF64}},ComplexF64,ComplexF64})
    Base.precompile(Tuple{typeof(mul!!),Matrix{ComplexF64},Transpose{ComplexF64, Matrix{ComplexF64}},Matrix{ComplexF64},ComplexF64,ComplexF64})
    Base.precompile(Tuple{typeof(mul!!),Matrix{ComplexF64},Transpose{ComplexF64, Matrix{ComplexF64}},Transpose{ComplexF64, Matrix{ComplexF64}},ComplexF64,ComplexF64})
    Base.precompile(Tuple{typeof(set_eltype),Type{Combiner},Type})
    Base.precompile(Tuple{typeof(set_eltype),Type{EmptyStorage{EmptyNumber, Dense{EmptyNumber, Vector{EmptyNumber}}}},Type})
    Base.precompile(Tuple{typeof(set_eltype),Type{IsWrappedArray{Combiner}},Type{Combiner},Type})
    Base.precompile(Tuple{typeof(tensor),Dense{ComplexF64, Vector{ComplexF64}},Vararg{Any}})
    Base.precompile(Tuple{typeof(tensor),Diag{Float64, Vector{Float64}},Vararg{Any}})
end
