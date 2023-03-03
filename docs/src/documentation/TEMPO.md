# Time-Evolving Matrix Product Operators

Time-Evolving Matrix Product Operators (TEMPO) uses matrix product states and operators to reduces the computational complexity and storage requirement for doing a path integral simulation. This allows for simulation of significantly larger systems with longer lengths of non-Markovian memory. The fundamental ideas are outlined in [TEMPO](https://dx.doi.org/10.1038/s41467-018-05617-3). The implementation closely follows the details presented in [TNPI](https://arxiv.org/abs/2106.12523).

## API
```@docs
TEMPO.propagate
```

```@docs
TEMPO.build_augmented_propagator
```