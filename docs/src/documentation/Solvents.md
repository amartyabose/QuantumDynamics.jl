# Solvents

For methods that simulate the solvent exactly using classical trajectories, we have a generic infrastructure developed that can be automatically leveraged.

## API
A generic solvent is encoded in the `Solvents` submodule:
```@docs
Solvents
```

The `Solvents` submodule defines a general `Solvent` type and an associated `PhaseSpace` type.
```@docs
Solvents.Solvent
```

```@docs
Solvents.PhaseSpace
```

Every solvent allows a phase-space point to be propagated under the internal forces that are present within the solvent and some external force (presumably coming from the system state) using
```@docs
Solvents.propagate_forced_bath
```

The internal forces of the solvent itself can be calculated using
```@docs
Solvents.bath_force
```
