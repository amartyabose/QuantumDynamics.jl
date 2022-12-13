# Quantum Dynamics

QuantumDynamics is an open-source software for the simulation of open quantum systems. Though written with performance in mind, QuantumDynamics provides a high throughput platform for experimentation with state-of-the-art approaches to method development.

The primary problem that QuantumDynamics is aimed at solving is simulation of the dynamics of a small quantum system coupled to a dissipative environment. Such a system-solvent decomposed problem can typically be represented by the Hamiltonian:
```math
\hat{H} = \hat{H}_0 + \hat{H}_\text{env}
```
where ``\hat{H}_0`` is the Hamiltonian of the isolated system and ``\hat{H}_\text{env}`` is the Hamiltonian corresponding to the environment and the interaction between the system and the environment.

## Installation
The QuantumDynamics package has not yet been registered. For the time being, the installation procedure directly uses the github repository. This can either be done by going into the package manager mode for Julia

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