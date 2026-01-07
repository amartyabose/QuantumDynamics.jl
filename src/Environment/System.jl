"General description of the system, particularly for trajectory based methods."
module Systems

using .Solvents: PhaseSpace

"Abstract type for all systems.
Every system needs to implement `Base.iterate` which either returns
the next sample point of the system, or the system _and_ the solvent."
abstract type System end

abstract type SystemPhaseSpace <: PhaseSpace end

"Abstract type for methods that propagate system and bath together."
abstract type CompositeSystem <: System end

end
