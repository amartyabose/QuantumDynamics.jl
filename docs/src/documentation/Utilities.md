# Utilities

Collection of some utilities for simulations.

## API

### ITensor Utilities

```@docs
Utilities.convert_ITensor_to_matrix
```

```@docs
Utilities.identity_MPO
```

```@docs
Utilities.MPO_to_MPS
```

```@docs
Utilities.MPS_to_MPO
```

```@docs
ITensors.expect
```

### Generic Utilities

```@docs
Utilities.get_BLAS_implementation
```

```@docs
Utilities.trapezoid
```

```@docs
Utilities.commutator
```

```@docs
Utilities.calculate_Liouvillian
```

```@docs
Utilities.ExternalField
```

```@docs
Utilities.hash_path
```

```@docs
Utilities.unhash_path
```

```@docs
Utilities.unhash_path_blips
```

```@docs
Utilities.apply_propagator
```

There are a few utilities for creating specific kinds of Hamiltonians.
```@docs
Utilities.create_tls_hamiltonian
```
```@docs
Utilities.create_nn_hamiltonian
```

Many of the algorithms require extra, method-specific arguments. These are implemented as subtypes of `Utilities.ExtraArgs`.
```@docs
Utilities.ExtraArgs
```

```@docs
Utilities.DiffEqArgs
```