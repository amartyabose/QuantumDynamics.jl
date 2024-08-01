# Utilities

Collection of some utilities for simulations.

## ITensor Utilities

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

## HDF5 Utilities

```@docs
Utilities.create_and_select_group
```

```@docs
Utilities.check_or_insert_value
```

```@docs
Utilities.merge_into
```

## Generic Utilities

```@docs
Utilities.get_BLAS_implementation
```

```@docs
Utilities.trapezoid
```

```@docs
Utilities.fourier_transform
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
Utilities.density_matrix_to_vector
```

```@docs
Utilities.density_matrix_vector_to_matrix
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

##  Extra arguments
Many of the algorithms require extra, method-specific arguments. These are implemented as subtypes of `Utilities.ExtraArgs`.
```@docs
Utilities.ExtraArgs
```

```@docs
Utilities.DiffEqArgs
```

```@docs
Utilities.TensorNetworkArgs
```