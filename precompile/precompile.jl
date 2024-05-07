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
    Base.precompile(Tuple{Type{QuantumDynamics.BlochRedfield.Params},Matrix{Float64},Array{Float64, 4},Vector{Float64},Nothing})
    Base.precompile(Tuple{Type{QuantumDynamics.EtaCoefficients.EtaCoeffs},ComplexF64,ComplexF64,Vector{ComplexF64},Vector{ComplexF64},Vector{ComplexF64}})
    Base.precompile(Tuple{Type{QuantumDynamics.QuAPI.Path{ComplexF64}},Vector{Int64},ComplexF64,Int64})
    Base.precompile(Tuple{Type{QuantumDynamics.QuAPI.Path{ComplexF64}},Vector{UInt64},ComplexF64,Int64})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{Hamiltonian::Matrix{ComplexF64}, Jw::Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff{Float64}}, β::Float64, ρ0::Matrix{ComplexF64}, dt::Float64, ntimes::Int64, sys_ops::Vector{Matrix{ComplexF64}}},typeof(QuantumDynamics.BlochRedfield.propagate)})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{Hamiltonian::Matrix{ComplexF64}, dt::Float64, ntimes::Int64},typeof(QuantumDynamics.Propagators.calculate_bare_propagators)})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{algorithm::String},Type{QuantumDynamics.Utilities.TensorNetworkArgs}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{fbU::Array{ComplexF64, 3}, Jw::Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff{Float64}}, β::Float64, dt::Float64, ntimes::Int64, svec::Matrix{Float64}},typeof(QuantumDynamics.Blip.build_augmented_propagator)})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{fbU::Array{ComplexF64, 3}, Jw::Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff{Float64}}, β::Float64, dt::Float64, ntimes::Int64, svec::Matrix{Float64}},typeof(QuantumDynamics.PCTNPI.build_augmented_propagator)})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{fbU::Array{ComplexF64, 3}, Jw::Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff{Float64}}, β::Float64, ρ0::Matrix{ComplexF64}, dt::Float64, ntimes::Int64, kmax::Int64, svec::Matrix{Float64}, extraargs::QuantumDynamics.Utilities.TensorNetworkArgs},typeof(QuantumDynamics.TEMPO.propagate)})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{fbU::Array{ComplexF64, 3}, Jw::Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff{Float64}}, β::Float64, ρ0::Matrix{ComplexF64}, dt::Float64, ntimes::Int64, kmax::Int64, svec::Matrix{Float64}},typeof(QuantumDynamics.QuAPI.propagate)})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{fbU::Array{ComplexF64, 3}, Jw::Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff{Float64}}, β::Float64, ρ0::Matrix{ComplexF64}, dt::Float64, ntimes::Int64, rmax::Int64, path_integral_routine::typeof(QuantumDynamics.Blip.build_augmented_propagator), extraargs::QuantumDynamics.Blip.BlipArgs},typeof(QuantumDynamics.TTM.propagate)})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{fbU::Array{ComplexF64, 3}, Jw::Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff{Float64}}, β::Float64, ρ0::Matrix{ComplexF64}, dt::Float64, ntimes::Int64, rmax::Int64, path_integral_routine::typeof(QuantumDynamics.Blip.build_augmented_propagator_parallel), extraargs::QuantumDynamics.Blip.BlipArgs},typeof(QuantumDynamics.TTM.propagate)})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{fbU::Array{ComplexF64, 3}, Jw::Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff{Float64}}, β::Float64, ρ0::Matrix{ComplexF64}, dt::Float64, ntimes::Int64, rmax::Int64, path_integral_routine::typeof(QuantumDynamics.QuAPI.build_augmented_propagator), extraargs::QuantumDynamics.QuAPI.QuAPIArgs},typeof(QuantumDynamics.TTM.propagate)})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{fbU::Array{ComplexF64, 3}, Jw::Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff{Float64}}, β::Float64, ρ0::Matrix{ComplexF64}, dt::Float64, ntimes::Int64, rmax::Int64, path_integral_routine::typeof(QuantumDynamics.QuAPI.build_augmented_propagator_parallel), extraargs::QuantumDynamics.QuAPI.QuAPIArgs},typeof(QuantumDynamics.TTM.propagate)})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{fbU::Array{ComplexF64, 3}, Jw::Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff{Float64}}, β::Float64, ρ0::Matrix{ComplexF64}, dt::Float64, ntimes::Int64, rmax::Int64, path_integral_routine::typeof(QuantumDynamics.TEMPO.build_augmented_propagator), extraargs::QuantumDynamics.Utilities.TensorNetworkArgs},typeof(QuantumDynamics.TTM.propagate)})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{tmpprops::Array{ComplexF64, 3}, path::Vector{UInt64}, group_Δs::Matrix{Float64}, sbar::Matrix{Float64}, η::Vector{QuantumDynamics.EtaCoefficients.EtaCoeffs}, propagator_type::String, nsteps::Int64, sdim2::Int64, val1::Vector{ComplexF64}, valend::Vector{ComplexF64}, valjkp::Matrix{ComplexF64}},typeof(QuantumDynamics.Blip.get_total_amplitude)})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{ξ::Float64, ωc::Float64},Type{QuantumDynamics.SpectralDensities.ExponentialCutoff}})
    Base.precompile(Tuple{typeof(Core.kwcall),@NamedTuple{ϵ::Float64, Δ::Float64},typeof(QuantumDynamics.Utilities.create_tls_hamiltonian)})
    Base.precompile(Tuple{typeof(QuantumDynamics.BlochRedfield.get_Rtensor),Vector{Float64},Matrix{ComplexF64},Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff{Float64}},Vector{Matrix{ComplexF64}},Float64})
    Base.precompile(Tuple{typeof(QuantumDynamics.SpectralDensities.reorganization_energy),QuantumDynamics.SpectralDensities.ExponentialCutoff{Float64}})
    Base.precompile(Tuple{typeof(QuantumDynamics.Utilities.commutator),Matrix{ComplexF64},Matrix{ComplexF64}})
    Base.precompile(Tuple{typeof(QuantumDynamics.Utilities.convert_ITensor_to_matrix),ITensor,Index{Int64},Index{Int64}})
end
