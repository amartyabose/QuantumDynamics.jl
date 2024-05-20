# Influence Functional Coefficients
The influence functional coefficients are discretizations of the bath correlation function required for simulations using the QuAPI influence functionals or QCPI with harmonic backreaction. QuantumDynamics provides facilities for generating these coefficients and storing them in a way that utilizes the limited time-translational symmetry that they demonstrate. These coefficients have been listed in many papers, the first being [QuAPI1](https://dx.doi.org/10.1063/1.469508).

## API
```@docs
EtaCoefficients.EtaCoeffs
```

```@docs
EtaCoefficients.calculate_η
```

```@docs
EtaCoefficients.ZetaCoeffs
```

```@docs
EtaCoefficients.calculate_ζ
```