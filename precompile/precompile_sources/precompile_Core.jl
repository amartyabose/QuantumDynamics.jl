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
    Base.precompile(Tuple{typeof(Core.Compiler.abstract_call_known),Core.Compiler.NativeInterpreter,Any,Core.Compiler.ArgInfo,Core.Compiler.StmtInfo,Core.Compiler.InferenceState,Int64})
    Base.precompile(Tuple{typeof(Core.Compiler.return_type),Core.Compiler.NativeInterpreter,DataType})
end
