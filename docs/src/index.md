# Quantum Dynamics

QuantumDynamics is an open-source software for the simulation of open quantum systems. Though written with performance in mind, QuantumDynamics provides a high throughput platform for experimentation with state-of-the-art approaches to method development.

The primary problem that QuantumDynamics is aimed at solving is simulation of the dynamics of a small quantum system coupled to a dissipative environment. Such a system-solvent decomposed problem can typically be represented by the Hamiltonian:
```math
\hat{H} = \hat{H}_0 + \hat{H}_\text{env}
```
where ``\hat{H}_0`` is the Hamiltonian of the isolated system and ``\hat{H}_\text{env}`` is the Hamiltonian corresponding to the environment and the interaction between the system and the environment.

Under Gaussian response theory, a molecular solvent can be mapped to a bath of harmonic oscillators bi-linearly interacting with the system:
```math
\hat{H}_\text{env} = \sum_j \frac{p_j^2}{2} + \frac{1}{2}\omega_j^2\left(x_j - \frac{c_j}{\omega_j^2}\hat{s}\right)^2
```
where ``\hat{s}`` is the system operator that couples with the bath modes. In such a harmonic mapping, the frequencies, ``\omega_j`` and the couplings ``c_j`` are characterized by the spectral density, which is related to the energy-gap autocorrelation function of the molecular solvent.
```math
J(\omega) = \frac{\pi}{2} \sum_j \frac{c_j^2}{\omega_j}\delta\left(\omega-\omega_j\right)
```

The simulations aim at computing the reduced density matrix corresponding to the system. The tracing out of the environment modes leads to temporally non-local interactions and non-Markovian dynamics. The non-Markovian interactions can be computed using the Feynman-Vernon influence functional, which when discretized give rise to the ``\eta``-coefficients used in this package.

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