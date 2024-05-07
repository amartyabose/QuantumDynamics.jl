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
    Base.precompile(Tuple{typeof(HDF5.API.__init__)})
    isdefined(HDF5, Symbol("#120#126")) && Base.precompile(Tuple{getfield(HDF5, Symbol("#120#126"))})
    isdefined(HDF5, Symbol("#121#127")) && Base.precompile(Tuple{getfield(HDF5, Symbol("#121#127"))})
    isdefined(HDF5, Symbol("#133#136")) && Base.precompile(Tuple{getfield(HDF5, Symbol("#133#136"))})
    isdefined(HDF5, Symbol("#134#137")) && Base.precompile(Tuple{getfield(HDF5, Symbol("#134#137"))})
end
