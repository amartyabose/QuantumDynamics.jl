# Quantum Dynamics

| **Documentation** | **Tests** |
|:-----------------:|:---------:|
|[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://amartyabose.github.io/QuantumDynamics.jl/dev/)|[![Run tests](https://github.com/amartyabose/QuantumDynamics.jl/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/amartyabose/QuantumDynamics.jl/actions/workflows/test.yml)|

QuantumDynamics is an open-source software for the simulation of open quantum systems. Though written with performance in mind, QuantumDynamics provides a high throughput platform for experimentation with state-of-the-art approaches to method development.

The primary problem that QuantumDynamics is aimed at solving is the simulation of the dynamics of a relatively small quantum system coupled to a dissipative environment. Such a system-solvent decomposed problem can typically be represented by the Hamiltonian:
```math
\hat{H} = \hat{H}_0 + \hat{H}_\text{env}
```
where ``\hat{H}_0`` is the Hamiltonian of the isolated system and ``\hat{H}_\text{env}`` is the Hamiltonian corresponding to the environment and the interaction between the system and the environment.

As demonstrated in the tutorials and the example codes, QuantumDynamics provides some approximate methods for simulating the dynamics of the system. However, the goal of this package is to provide access to more state-of-the-art techniques based on path integrals, tensor networks and other ideas in such a manner that all of these methods can be used as far as possible in a composable manner.

## Installation
The QuantumDynamics.jl package is registered. The installation can either be done by going into the Pkg REPL mode for Julia

```bash
~ julia
```

```
julia> ]
pkg> add QuantumDynamics
```

or by using the `Pkg` package manager in a script:

```julia
julia> using Pkg
julia> Pkg.add("QuantumDynamics")
```

## Citation
If you use QuantumDynamics in your work, please cite the [QuantumDynamics.jl paper](https://pubs.aip.org/aip/jcp/article/158/20/204113/2892511/QuantumDynamics-jl-A-modular-approach-to):
```bibtex
@article{10.1063/5.0151483,
    author = {Bose, Amartya},
    title = "{QuantumDynamics.jl: A modular approach to simulations of dynamics of open quantum systems}",
    journal = {The Journal of Chemical Physics},
    volume = {158},
    number = {20},
    year = {2023},
    month = {05},
    abstract = "{A simulation of the non-adiabatic dynamics of a quantum system coupled to dissipative environments poses significant challenges. New sophisticated methods are regularly being developed with an eye toward moving to larger systems and more complicated descriptions of solvents. Many of these methods, however, are quite difficult to implement and debug. Furthermore, trying to make the individual algorithms work together through a modular application programming interface can be quite difficult as well. We present a new, open-source software framework, QuantumDynamics.jl, designed to address these challenges. It provides implementations of a variety of perturbative and non-perturbative methods for simulating the dynamics of these systems. Most prominently, QuantumDynamics.jl supports hierarchical equations of motion and methods based on path integrals. An effort has been made to ensure maximum compatibility of the interface between the various methods. Additionally, QuantumDynamics.jl, being built on a high-level programming language, brings a host of modern features to explorations of systems, such as the usage of Jupyter notebooks and high level plotting, the possibility of leveraging high-performance machine learning libraries for further development. Thus, while the built-in methods can be used as end-points in themselves, the package provides an integrated platform for experimentation, exploration, and method development.}",
    issn = {0021-9606},
    doi = {10.1063/5.0151483},
    url = {https://doi.org/10.1063/5.0151483},
    note = {204113},
    eprint = {https://pubs.aip.org/aip/jcp/article-pdf/doi/10.1063/5.0151483/17794821/204113\_1\_5.0151483.pdf},
}
```