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
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:Hamiltonian, :Jw, :β, :ρ0, :dt, :ntimes, :sys_ops), Tuple{Matrix{ComplexF64}, Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff}, Float64, Matrix{ComplexF64}, Float64, Int64, Vector{Matrix{ComplexF64}}}},typeof(QuantumDynamics.BlochRedfield.propagate)})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:Hamiltonian, :dt, :ntimes), Tuple{Matrix{ComplexF64}, Float64, Int64}},typeof(QuantumDynamics.Propagators.calculate_bare_propagators)})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:fbU, :Jw, :β, :dt, :ntimes, :svec), Tuple{Array{ComplexF64, 3}, Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff}, Float64, Float64, Int64, Matrix{Float64}}},typeof(QuantumDynamics.Blip.build_augmented_propagator)})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:fbU, :Jw, :β, :dt, :ntimes, :svec), Tuple{Array{ComplexF64, 3}, Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff}, Float64, Float64, Int64, Matrix{Float64}}},typeof(QuantumDynamics.PCTNPI.build_augmented_propagator)})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:fbU, :Jw, :β, :ρ0, :dt, :ntimes, :kmax, :svec), Tuple{Array{ComplexF64, 3}, Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff}, Float64, Matrix{ComplexF64}, Float64, Int64, Int64, Matrix{Float64}}},typeof(QuantumDynamics.QuAPI.propagate)})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:fbU, :Jw, :β, :ρ0, :dt, :ntimes, :kmax, :svec), Tuple{Array{ComplexF64, 3}, Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff}, Float64, Matrix{ComplexF64}, Float64, Int64, Int64, Matrix{Float64}}},typeof(QuantumDynamics.TEMPO.propagate)})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:fbU, :Jw, :β, :ρ0, :dt, :ntimes, :rmax, :path_integral_routine, :extraargs), Tuple{Array{ComplexF64, 3}, Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff}, Float64, Matrix{ComplexF64}, Float64, Int64, Int64, typeof(QuantumDynamics.TEMPO.build_augmented_propagator), QuantumDynamics.TEMPO.TEMPOArgs}},typeof(QuantumDynamics.TTM.propagate)})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:ξ, :ωc), Tuple{Float64, Float64}},Type{QuantumDynamics.SpectralDensities.ExponentialCutoff}})
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:ϵ, :Δ), Tuple{Float64, Float64}},typeof(QuantumDynamics.Utilities.create_tls_hamiltonian)})
    Base.precompile(Tuple{typeof(QuantumDynamics.BlochRedfield.get_Rtensor),Vector{Float64},Matrix{ComplexF64},Vector{QuantumDynamics.SpectralDensities.ExponentialCutoff},Vector{Matrix{ComplexF64}},Float64})
    Base.precompile(Tuple{typeof(QuantumDynamics.EtaCoefficients.calculate_η),Vector{Float64},Vector{Float64},Float64,Float64,Int64,Bool,Bool,Bool})
    Base.precompile(Tuple{typeof(QuantumDynamics.Utilities.commutator),Matrix{ComplexF64},Matrix{ComplexF64}})
    Base.precompile(Tuple{typeof(QuantumDynamics.Utilities.convert_ITensor_to_matrix),ITensor,Index{Int64},Index{Int64}})
end
