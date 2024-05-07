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
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{current_time::ComplexF64, nsite::Int64, reverse_step::Bool, sweep::Int64, observer!::NoObserver, normalize::Bool, maxdim::Int64, mindim::Int64, cutoff::Float64, noise::Int64},typeof(sub_sweep_update),Base.Order.ReverseOrdering{Base.Order.ForwardOrdering},Function,ProjMPOApply,ComplexF64,MPS})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{current_time::ComplexF64, outputlevel::Int64, time_step::ComplexF64, normalize::Bool, direction::Base.Order.ForwardOrdering, noise::Int64, which_decomp::Nothing, svd_alg::Nothing, cutoff::Float64, maxdim::Int64, mindim::Int64, maxtruncerr::Float64},typeof(region_update!),Val{2},Val{false},Function,ProjMPOApply,MPS,Int64})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{current_time::ComplexF64, outputlevel::Int64, time_step::ComplexF64, normalize::Bool, direction::Base.Order.ReverseOrdering{Base.Order.ForwardOrdering}, noise::Int64, which_decomp::Nothing, svd_alg::Nothing, cutoff::Float64, maxdim::Int64, mindim::Int64, maxtruncerr::Float64},typeof(region_update!),Val{2},Val{false},Function,ProjMPOApply,MPS,Int64})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{current_time::Int64, nsite::Int64, reverse_step::Bool, sweep::Int64, observer!::NoObserver, normalize::Bool, maxdim::Int64, mindim::Int64, cutoff::Float64, noise::Int64},typeof(sub_sweep_update),Base.Order.ForwardOrdering,Function,ProjMPOApply,ComplexF64,MPS})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{current_time::Int64, outputlevel::Int64, time_step::ComplexF64, normalize::Bool, direction::Base.Order.ForwardOrdering, noise::Int64, which_decomp::Nothing, svd_alg::Nothing, cutoff::Float64, maxdim::Int64, mindim::Int64, maxtruncerr::Float64},typeof(region_update!),Val{2},Val{false},Function,ProjMPOApply,MPS,Int64})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{nsite::Int64, current_time::Int64, reverse_step::Bool, sweep::Int64, observer!::NoObserver, normalize::Bool, maxdim::Int64, mindim::Int64, cutoff::Float64, noise::Int64},typeof(sweep_update),TDVPOrder{2, Base.Order.ForwardOrdering()},Function,ProjMPOApply,ComplexF64,MPS})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{nsite::Int64, reverse_step::Bool, current_time::ComplexF64, outputlevel::Int64, time_step::ComplexF64, normalize::Bool, direction::Base.Order.ForwardOrdering, noise::Int64, which_decomp::Nothing, svd_alg::Nothing, cutoff::Float64, maxdim::Int64, mindim::Int64, maxtruncerr::Float64},typeof(region_update!),Function,ProjMPOApply,MPS,Int64})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{nsweeps::Int64, reverse_step::Bool, cutoff::Float64, maxdim::Int64},typeof(alternating_update),Function,ProjMPOApply,ComplexF64,MPS})
end
