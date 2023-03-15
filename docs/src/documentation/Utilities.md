# Utilities

Collection of some utilities for simulations.

## API
```@docs
Utilities.get_BLAS_implementation
```

```@docs
Utilities.ExternalField
```

```@docs
Utilities.unhash_path
```

```@docs
Utilities.apply_propagator
```

```@docs
Utilities.convert_ITensor_to_matrix
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