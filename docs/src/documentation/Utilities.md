# Utilities

Collection of some utilities for simulations.

## API
```@docs
Utilities.ExternalField
```

```@docs
Utilities.unhash_path
```

```@docs
Utilities.apply_propagator
```

Many of the algorithms require extra, method-specific arguments. These are implemented as subtypes of `Utilities.ExtraArgs`.
```@docs
Utilities.ExtraArgs
```

```@docs
Utilities.DiffEqArgs
```