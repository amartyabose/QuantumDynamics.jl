function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    # QuAPI
    Base.precompile(Tuple{Type{QuAPI.Path{ComplexF64}},Vector{Int64},Any,Int64})   # time: 0.002645703
    Base.precompile(Tuple{Type{QuAPI.QuAPIArgs}})   # time: 0.001081874

    # Blip
    isdefined(QuantumDynamics.Blip, Symbol("#7#17")) && Base.precompile(Tuple{getfield(QuantumDynamics.Blip, Symbol("#7#17")),Int64})   # time: 0.05048021
    Base.precompile(Tuple{typeof(Blip.setup_simulation),Matrix{Float64}})   # time: 0.013447023
    Base.precompile(Tuple{typeof(Core.kwcall),Core.NamedTuple{(:tmpprops, :path, :group_Δs, :sbar, :η, :propagator_type, :nsteps, :sdim2, :val1, :valend, :valjkp), <:Tuple{Core.Array{Base.ComplexF64, 3}, Core.Any, Base.Matrix{Core.Float64}, Base.Matrix{Core.Float64}, Base.Vector{QuantumDynamics.EtaCoefficients.EtaCoeffs}, Core.String, Core.Int64, Core.Int64, Base.Vector{Base.ComplexF64}, Base.Vector{Base.ComplexF64}, Base.Matrix{Base.ComplexF64}}},typeof(Blip.get_total_amplitude)})   # time: 0.007939419
    Base.precompile(Tuple{Type{Blip.BlipArgs}})   # time: 0.001180291

    # BlochRedfield
    Base.precompile(Tuple{typeof(BlochRedfield.func_BRME),Matrix{ComplexF64},BlochRedfield.Params,Float64})   # time: 0.11660305
    Base.precompile(Tuple{typeof(BlochRedfield.func_BRME),Any,BlochRedfield.Params,Float64})   # time: 0.03534429
    Base.precompile(Tuple{typeof(BlochRedfield.func_BRME),Matrix{ComplexF64},Any,Any})   # time: 0.003835044

    # EtaCoefficients
    Base.precompile(Tuple{typeof(EtaCoefficients.calculate_η),Vector{Float64},Vector{Float64},Float64,Float64,Int64,Bool,Bool,Bool})   # time: 0.42108637

    # Propagators
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{Hamiltonian::Matrix{ComplexF64}, dt::Float64, ntimes::Int64},typeof(Propagators.calculate_bare_propagators)})   # time: 0.4211968

    # SpectralDensities
    Base.precompile(Tuple{typeof(SpectralDensities.reorganization_energy),SpectralDensities.ExponentialCutoff})   # time: 0.3470814
    Base.precompile(Tuple{typeof(SpectralDensities.eval_spectrum),SpectralDensities.ExponentialCutoff,Real,Float64})   # time: 0.003393833

    # TEMPO
    isdefined(QuantumDynamics.TEMPO, Symbol("#2#9")) && Base.precompile(Tuple{getfield(QuantumDynamics.TEMPO, Symbol("#2#9")),Int64})   # time: 0.001370708

    # TTM
    Base.precompile(Tuple{typeof(TTM.get_Ts),Array{<:Complex, 3}})   # time: 0.02414457
    Base.precompile(Tuple{typeof(TTM.get_Ts),Array{ComplexF64, 3}})   # time: 0.005978883

    # Utilities
    Base.precompile(Tuple{typeof(Utilities.build_path_amplitude_mps),Matrix{ComplexF64},Vector{Index{Int64}}})   # time: 0.07074379
    Base.precompile(Tuple{typeof(Utilities.apply_contract_propagator),MPS,MPO})   # time: 0.028912203
    Base.precompile(Tuple{typeof(Utilities.unhash_path),Int64,Vector{UInt64},Int64})   # time: 0.018552251
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{propagators::Array{ComplexF64, 3}, ρ0::Matrix{ComplexF64}, ntimes::Int64, dt::Float64},typeof(Utilities.apply_propagator)})   # time: 0.015000333
    Base.precompile(Tuple{typeof(Utilities.convert_ITensor_to_matrix),Any,Index{Int64},Index{Int64}})   # time: 0.012299342
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{ϵ::Float64, Δ::Float64},typeof(Utilities.create_tls_hamiltonian)})   # time: 0.011816128
    Base.precompile(Tuple{typeof(Utilities.extend_path_amplitude_mps),MPS,Matrix{ComplexF64},Vector{Index{Int64}}})   # time: 0.007114587
    Base.precompile(Tuple{typeof(Utilities.extend_path_amplitude_mps_beyond_memory),MPS,Matrix{ComplexF64},Vector{Index{Int64}}})   # time: 0.004667583
    Base.precompile(Tuple{typeof(Utilities.convert_ITensor_to_matrix),ITensor,Index{Int64},Index{Int64}})   # time: 0.004362336
    Base.precompile(Tuple{Type{Utilities.TensorNetworkArgs}})   # time: 0.003767875
    Base.precompile(Tuple{typeof(Utilities.commutator),Matrix{ComplexF64},Matrix{ComplexF64}})   # time: 0.002387746
    Base.precompile(Tuple{typeof(Utilities.blip_dist_criterion),Vector{UInt64},Int64})   # time: 0.002369574
    Base.precompile(Tuple{Type{Utilities.DiffEqArgs}})   # time: 0.001335293
    Base.precompile(Tuple{typeof(Utilities.blip_dist_criterion),Any,Int64})   # time: 0.001046168

    # BMatrix
    Base.precompile(Tuple{typeof(BMatrix.get_B_matrix),Vector{Float64},Vector{Float64},Float64,Float64,Int64})   # time: 0.05446195

    # ComplexPISetup
    Base.precompile(Tuple{typeof(ComplexPISetup.get_complex_time_propagator),Matrix{ComplexF64},Float64,Float64,Int64})   # time: 0.006184968

    # ComplexTNPI
    Base.precompile(Tuple{typeof(ComplexTNPI.get_Bmat_MPO_left),Matrix{Float64},Vector,Int64,Vector{Matrix{ComplexF64}},Int64})   # time: 0.017964996
    Base.precompile(Tuple{typeof(ComplexTNPI.get_Bmat_MPO_right),Matrix{Float64},Vector{Index{Int64}},Int64,Vector{Matrix{ComplexF64}},Int64})   # time: 0.015762536
    Base.precompile(Tuple{typeof(ComplexTNPI.get_Bmat_MPO_left),Matrix{Float64},Vector{Index{Int64}},Int64,Vector{Matrix{ComplexF64}},Int64})   # time: 0.015465842
    Base.precompile(Tuple{typeof(ComplexTNPI.get_Bmat_MPO_right),Matrix{Float64},Vector,Int64,Vector{Matrix{ComplexF64}},Int64})   # time: 0.008204346
end
