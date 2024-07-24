function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    isdefined(QuantumDynamics.Blip, Symbol("#7#17")) && Base.precompile(Tuple{getfield(QuantumDynamics.Blip, Symbol("#7#17")),Int64})   # time: 0.053283364
    Base.precompile(Tuple{typeof(Core.kwcall),Core.NamedTuple{(:tmpprops, :path, :group_Δs, :sbar, :η, :propagator_type, :nsteps, :sdim2, :val1, :valend, :valjkp), <:Tuple{Core.Array{Base.ComplexF64, 3}, Core.Any, Base.Matrix{Core.Float64}, Base.Matrix{Core.Float64}, Base.Vector{QuantumDynamics.EtaCoefficients.EtaCoeffs}, Core.String, Core.Int64, Core.Int64, Base.Vector{Base.ComplexF64}, Base.Vector{Base.ComplexF64}, Base.Matrix{Base.ComplexF64}}},typeof(get_total_amplitude)})   # time: 0.008192462
end
