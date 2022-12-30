# η-Coefficients
The η-coefficients [makriTensorPropagatorIterativeI1995](@cite) are discretizations of the bath correlation function required for simulations using the QuAPI influence functionals. QuantumDynamics provides facilities for generating these coefficients and storing them in a way that utilizes the limited time-translational symmetry that they demonstrate. 

## API
```@docs
EtaCoefficients.EtaCoeffs
```

```@docs
EtaCoefficients.calculate_η
```