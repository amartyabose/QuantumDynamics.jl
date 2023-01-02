# Quantum Dynamics

| **Documentation** |
|:-----------------:|
|[![docs-dev][docsdev-img]][docsdev-url]|

[docsdev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docsdev-url]: https://amartyabose.github.io/QuantumDynamics/dev/

QuantumDynamics is an open-source software for the simulation of open quantum systems. Though written with performance in mind, QuantumDynamics provides a high throughput platform for experimentation with state-of-the-art approaches to method development.

The primary problem that QuantumDynamics is aimed at solving is simulation of the dynamics of a relatively small quantum system coupled to a dissipative environment. Such a system-solvent decomposed problem can typically be represented by the Hamiltonian:
```math
\hat{H} = \hat{H}_0 + \hat{H}_\text{env}
```
where ``\hat{H}_0`` is the Hamiltonian of the isolated system and ``\hat{H}_\text{env}`` is the Hamiltonian corresponding to the environment and the interaction between the system and the environment.

As demonstrated in the tutorials and the example codes, QuantumDynamics provides some approximate methods for simulating the dynamics of the system. However, the goal of this package is to provide access to more state-of-the-art techniques based on path integrals, tensor networks and other ideas in such a manner that all of these methods can be used as far as possible in a composable manner.

## Installation
The QuantumDynamics package has not yet been registered. For the time being, the installation procedure directly uses the github repository. This can either be done by going into the Pkg REPL mode for Julia

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