# Quantum Dynamics

QuantumDynamics is an open-source software for the simulation of open quantum systems. Though written with performance in mind, QuantumDynamics provides a high throughput platform for experimentation with state-of-the-art approaches to method development.

## Installation
The QuantumDynamics package has not yet been registered. For the time being, the installation procedure directly uses the github repository. This can either be done by going into the package manager mode for julia

```bash
~ julia
```

```
julia> ]
pkg> add https://github.com/amartyabose/QuantumDynamics
```

or by using the `Pkg` package manager in a script:

```julia
julia> using Pkg
julia> Pkg.add("https://github.com/amartyabose/QuantumDynamics")
```