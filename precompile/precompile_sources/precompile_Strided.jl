function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(Base.Broadcast.dotview),StridedView{<:Any, 3},Union{Colon, AbstractRange{<:Union{Int128, Int16, Int32, Int64, Int8, UInt128, UInt16, UInt32, UInt64, UInt8}}, Int128, Int16, Int32, Int64, Int8, UInt128, UInt16, UInt32, UInt64, UInt8},Colon,Colon})   # time: 0.008667526
end
