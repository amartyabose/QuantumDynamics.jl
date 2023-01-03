# Quasi-Adiabatic Propagator Path Integral (QuAPI)

This module provides the basic interface for simulating a system using QuAPI. Though the implementation does not follow the algorithm in the original papers, the first papers to outline the method are [MakriQuAPI1](https://dx.doi.org/10.1063/1.469508) and [MakriQuAPI2](https://dx.doi.org/10.1063/1.469509). For an overall review of the ideas involved, consider reading [MakriNumericalPathIntegrals](https://doi.org/10.1063/1.531046). 

## API
The API has three important end-points for user interface.

First, there is the `propagate` function for propagating a given reduced density matrix.

```@docs
QuAPI.propagate
```

Then, there is the `build_augmented_propagator` function for computing the augmented propagator incorporating the solvent effects through the Feynman-Vernon influence functional. This currently only does a full path calculation and does not iterate.

```@docs
QuAPI.build_augmented_propagator
```

QuAPI allows for path filtering based on the absolute value of the amplitude of a path. This cutoff threshold is specified using the `QuAPIArgs` structure. `QuAPI.propagate` and `QuAPI.build_augmented_propagator` use objects of this structure, with a default value being the default constructed objected.

```@docs
QuAPI.QuAPIArgs
```