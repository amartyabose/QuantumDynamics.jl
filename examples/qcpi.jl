using Revise
using QuantumDynamics
using LaTeXStrings
import PyPlot;
const plt = PyPlot;
using Random

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254

function new_figure(plot_type="full")
    if plot_type == "full"
        fig = plt.figure(; figsize=(3.175, 2.25))
    elseif plot_type == "half"
        fig = plt.figure(; figsize=(1.56, 1.4))
    end
    ax = fig.add_subplot(111)

    for axis in ["top", "bottom", "left", "right"]
        ax.spines[axis].set_linewidth(0.1)
        ax.spines[axis].set_color("gray")
    end
    ax.tick_params(width=0.1, direction="in", color="gray")

    fig, ax
end

function qcpi(num_points)
    # specify the system Hamiltonian
    H0 = Matrix{ComplexF64}([
        1.0 -1.0
        -1.0 -1.0
    ])

    # specify the spectral density and the inverse temperature
    Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)
    β = 5.0

    ρ0 = Matrix{ComplexF64}([
        1.0 0.0
        0.0 0.0
    ])

    dt = 0.25
    ntimes = 100

    # obtain the basic TTM results
    @info "Starting TTM calculation"
    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt, ntimes)
    @time times, ρs = TTM.propagate(; fbU=fbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=9, extraargs=TEMPO.TEMPOArgs(), path_integral_routine=TEMPO.build_augmented_propagator)

    # discretize the spectral density and create a harmonic bath solvent
    # for an atomistic solvent, here we would use the actual description based on an appropriate force field or ab initio DFT calculation
    ω, c = SpectralDensities.discretize(Jw, 100)
    hb = Solvents.HarmonicBath(β, [ω], [c], [1.0 -1.0], num_points)

    # η = EtaCoefficients.calculate_η(Jw; β, dt, kmax=5, imaginary_only=true)
    # ζ = EtaCoefficients.calculate_ζ(Jw; dt, kmax=5)
    # display([-2*imag.(η.ηmn) ζ.ζmn[:,1].+ζ.ζmn[:,2]])
    # @assert -2*imag.(η.ηmn) ≈ ζ.ζmn[:,1].+ζ.ζmn[:,2]

    # display([-2*imag.(η.ηmm)+dt*SpectralDensities.reorganization_energy(Jw) ζ.ζmm[1].+ζ.ζmm[2]])
    # display([-2*imag.(η.η00)+0.25*dt*SpectralDensities.reorganization_energy(Jw) ζ.ζ00])

    # calculate EACP dynamics
    # @info "Starting fixed-EACP calculation"
    # @time EACP_fbU = Propagators.calculate_average_reference_propagators(; Hamiltonian=H0, solvent=hb, classical_dt=dt / 100, dt, ntimes, ρ0, reference_choice="fixed", verbose=true)
    # @info "Reference propagators calculated"
    # @time times_EACP, ρs_EACP = Utilities.apply_propagator(; propagators=EACP_fbU, ρ0, ntimes, dt)
    # @info "fixed-EACP done"

    # @info "Starting Ehrenfest calculation"
    # @time EACP_fbU = Propagators.calculate_average_reference_propagators(; Hamiltonian=H0, solvent=hb, classical_dt=dt / 100, dt, ntimes, ρ0, reference_choice="ehrenfest", verbose=true)
    # @info "Reference propagators calculated"
    # @time times_ehrenfest, ρs_ehrenfest = Utilities.apply_propagator(; propagators=EACP_fbU, ρ0, ntimes, dt)
    # @info "Ehrenfest done"

    # @info "Starting DCSH calculation"
    # @time EACP_fbU = Propagators.calculate_average_reference_propagators(; Hamiltonian=H0, solvent=hb, classical_dt=dt / 100, dt, ntimes, ρ0, reference_choice="dcsh", verbose=true)
    # @info "Reference propagators calculated"
    # @time times_dcsh, ρs_dcsh = Utilities.apply_propagator(; propagators=EACP_fbU, ρ0, ntimes, dt)
    # @info "DCSH done"

    # @info "Starting max_prob calculation"
    # @time EACP_fbU = Propagators.calculate_average_reference_propagators(; Hamiltonian=H0, solvent=hb, classical_dt=dt / 100, dt, ntimes, ρ0, reference_choice="max_prob", verbose=true)
    # @info "Reference propagators calculated"
    # @time times_max_prob, ρs_max_prob = Utilities.apply_propagator(; propagators=EACP_fbU, ρ0, ntimes, dt)
    # @info "max_prob done"

    # simulate QCPI

    @info "Starting QCPI"
    # Random.seed!(261188)
    @time times_QCPI_fixed, ρs_QCPI_fixed = QCPI.propagate(; Hamiltonian=H0, Jw=[Jw], solvent=hb, ρ0, classical_dt=dt / 100, dt, ntimes, kmax=5, reference_choice="fixed", extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate, verbose=true)
    @info "QCPI done"
    @info "Starting QCPI DCSH"
    # Random.seed!(261188)
    @time times_QCPI_dcsh, ρs_QCPI_dcsh = QCPI.propagate(; Hamiltonian=H0, Jw=[Jw], solvent=hb, ρ0, classical_dt=dt / 100, dt, ntimes, kmax=5, reference_choice="dcsh", extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate_td_reference, verbose=true)
    @info "QCPI done"
    @info "Starting QCPI ehrenfest"
    # Random.seed!(261188)
    @time times_QCPI_ehrenfest, ρs_QCPI_ehrenfest = QCPI.propagate(; Hamiltonian=H0, Jw=[Jw], solvent=hb, ρ0, classical_dt=dt / 100, dt, ntimes, kmax=5, reference_choice="ehrenfest", extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate_td_reference, verbose=true)
    @info "QCPI done"
    @info "Starting QCPI max_prob"
    # Random.seed!(261188)
    @time times_QCPI_max_prob, ρs_QCPI_max_prob = QCPI.propagate(; Hamiltonian=H0, Jw=[Jw], solvent=hb, ρ0, classical_dt=dt / 100, dt, ntimes, kmax=5, reference_choice="max_prob", extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate_td_reference, verbose=true)
    @info "QCPI done"

    new_figure()
    plt.plot(times, real.(ρs[:, 1, 1]), lw=0.75, label="Basic Path Integral")
    plt.plot(times_QCPI_fixed, real.(ρs_QCPI_fixed[:, 1, 1]), lw=0.75, label="Fixed")
    plt.plot(times_QCPI_dcsh, real.(ρs_QCPI_dcsh[:, 1, 1]), lw=0.75, label="DCSH")
    plt.plot(times_QCPI_ehrenfest, real.(ρs_QCPI_ehrenfest[:, 1, 1]), lw=0.75, label="Ehrenfest")
    plt.plot(times_QCPI_max_prob, real.(ρs_QCPI_max_prob[:, 1, 1]), lw=0.75, label="max_prob")
    plt.legend()
    plt.xlabel(L"\Omega t")
    plt.ylabel(L"P(t)")
    plt.savefig("qcpi_$(num_points).pdf"; bbox_inches="tight")
end

function qcpi_drude(num_points)
    # specify the system Hamiltonian
    H0 = Matrix{ComplexF64}([
        1.0 -1.0
        -1.0 -1.0
    ])

    # specify the spectral density and the inverse temperature
    Jw = SpectralDensities.DrudeLorentz(; λ=1.5, γ=7.5)
    β = 5.0

    ρ0 = Matrix{ComplexF64}([
        1.0 0.0
        0.0 0.0
    ])

    dt = 0.25
    ntimes = 100

    # obtain the basic TTM results
    @info "Starting TTM calculation"
    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt, ntimes)
    @time times, ρs = TTM.propagate(; fbU=fbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=9, extraargs=TEMPO.TEMPOArgs(), path_integral_routine=TEMPO.build_augmented_propagator)

    # discretize the spectral density and create a harmonic bath solvent
    # for an atomistic solvent, here we would use the actual description based on an appropriate force field or ab initio DFT calculation
    ω, c = SpectralDensities.discretize(Jw, 100)
    # display([ω c])
    hb = Solvents.HarmonicBath(β, [ω], [c], [1.0 -1.0], num_points)

    # calculate EACP dynamics
    @info "Starting EACP calculation"
    @time EACP_fbU = Propagators.calculate_average_reference_propagators(; Hamiltonian=H0, solvent=hb, classical_dt=dt / 100, dt, ntimes, ρ0)
    @info "Reference propagators calculated"
    @time times_EACP, ρs_EACP = Utilities.apply_propagator(; propagators=EACP_fbU, ρ0, ntimes, dt)
    @info "EACP done"

    # simulate QCPI
    @info "Starting QCPI"
    @time times_QCPI, ρs_QCPI = QCPI.propagate(; Hamiltonian=H0, Jw, solvent=hb, ρ0, classical_dt=dt / 100, dt, ntimes, kmax=5, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate)
    @info "QCPI done"

    new_figure()
    plt.plot(times, real.(ρs[:, 1, 1]), lw=0.75, label="Basic Path Integral")
    plt.plot(times_EACP, real.(ρs_EACP[:, 1, 1]), lw=0.75, label="EACP")
    plt.plot(times_QCPI, real.(ρs_QCPI[:, 1, 1]), lw=0.75, label="QCPI")
    plt.legend()
    plt.xlabel(L"\Omega t")
    plt.ylabel(L"P(t)")
    plt.savefig("qcpi_drude_$(num_points).pdf"; bbox_inches="tight")
end