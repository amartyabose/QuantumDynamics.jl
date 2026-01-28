"General description of the system, particularly for trajectory based methods."
module Systems

using ..Solvents: PhaseSpace
using ..SolventsX: Solvent
using LinearAlgebra: diag

abstract type SystemPhaseSpace <: PhaseSpace end

"""Abstract type for methods that propagate system and bath together.
All such systems must implement `Base.iterate` which returns the next
sample point of the system and the solvent."""
abstract type CompositeSystem end

"""Abstract type for methods that mapping Hamiltonian methods."""
abstract type MappedSystem <: CompositeSystem end

### Spin-mapping methods.

"Abstract type for all Stratonovich–Weyl transforms."
abstract type SWTransform end

abstract type WTransform <: SWTransform end
abstract type QTransform <: SWTransform end
abstract type PTransform <: SWTransform end

R²(::Type{WTransform}, d::Integer) = 2sqrt(d+1)
R²(::Type{QTransform}, d::Integer) = 2.0
R²(::Type{PTransform}, d::Integer) = 2.0 * (d+1)

γ(t::Type{<:SWTransform}, d::Integer) = (R²(t, d) - 2.0) / d
γ(::Type{QTransform}, d::Integer) = 0.0

dual(W::Type{WTransform}) = W
dual(::Type{PTransform}) = QTransform
dual(::Type{QTransform}) = PTransform

rescale_factor(s::Type{<:SWTransform}, s̄::Type{<:SWTransform}, d::Integer) = sqrt(R²(s̄, d) / R²(s, d))
rescale_factor(::Type{WTransform}, ::Type{WTransform}, ::Integer) = 1

"""Abstract type for all spin-mapped systems.
It should contain the following elements always:
  - `transform`: the Stratonovich–Weyl transform
  - `d`: the dimensionality of the system
  - `γₛ`: the zero-point energy parameter of the Stratonovich–Weyl transform
"""
abstract type SpinMappedSystem <: MappedSystem end

# TODO: Find a better place for this generic method.
"Perform the relevant SW transform of the operator `op`."
function transform_op(sys::SpinMappedSystem,
                      op::AbstractMatrix,
                      X::Vector{<:Real}, P::Vector{<:Real})
    opt = 0.5 * sum(diag(op) .* (X.^2 .+ P.^2 .- sys.γₛ))

    for n = 1:sys.d, m = n+1:sys.d
        opt += 0.5 * op[n,m] * (X[n] - im * P[n]) * (X[m] + im * P[m])
        opt += 0.5 * op[m,n] * (X[m] - im * P[m]) * (X[n] + im * P[n])
    end

    opt
end

end
