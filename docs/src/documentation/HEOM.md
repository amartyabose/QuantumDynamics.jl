# Hierarchy Equations of Motion

This module provides the necessary tools for doing HEOM simulations. While the equations of motion have been reported in many articles, the current implementation is based on a relatively recent paper, [HEOM1](https://doi.org/10.1021/acs.jctc.5b00488). It should be noted that this is a particularly naive implementation of the algorithm that only works for the simplest of case of a Drude-Lorentz spectral density.

## API
```@docs
HEOM.propagate
```