"General description of the system, particularly for trajectory based methods."
module Systems

using ..Solvents: PhaseSpace

abstract type SystemPhaseSpace <: PhaseSpace end

"""Abstract type for methods that propagate system and bath together.
All such systems must implement `Base.iterate` which returns the next
sample point of the system and the solvent."""
abstract type CompositeSystem end

end
