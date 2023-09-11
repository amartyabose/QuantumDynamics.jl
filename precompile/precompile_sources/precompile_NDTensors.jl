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
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{2}},Tuple{Int64, Int64, Int64},Tuple{Int64, Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{2}},Tuple{Int64, Int64},Tuple{Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{3}},NTuple{4, Int64},Tuple{Int64, Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{3}},NTuple{4, Int64},Tuple{Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{3}},Tuple{Int64, Int64, Int64},NTuple{4, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{3}},Tuple{Int64, Int64, Int64},Tuple{Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{3}},Tuple{Int64, Int64},Tuple{Int64, Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{4}},NTuple{5, Int64},Tuple{Int64, Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{4}},NTuple{5, Int64},Tuple{Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{4}},Tuple{Int64, Int64, Int64},NTuple{5, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{4}},Tuple{Int64, Int64, Int64},Tuple{Int64, Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{5}},NTuple{4, Int64},Tuple{Int64, Int64, Int64}})
    Base.precompile(Tuple{typeof(contract_labels),Type{Val{5}},Tuple{Int64, Int64, Int64},NTuple{4, Int64}})
end
