# Quantum-Classical Path Integral

QCPI provides a rigorous way of coupling a classical-like solvent to a quantum system. The method has been outlined in [QCPI1](https://dx.doi.org/10.1063/1.4767931), [QCPI2](https://dx.doi.org/10.1063/1.4767980), [reference propagators](https://dx.doi.org/10.1063/1.4767980).

## API
```@docs
QCPI.propagate
```

The solvent shown here is encoded in the `Solvents` submodule:
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