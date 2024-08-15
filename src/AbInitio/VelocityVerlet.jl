module VelocityVerlet

using LinearAlgebra

using AtomsIO
using Unitful, UnitfulAtomic

using ..Runners, ..MolecularUtilities

function velocity_verlet(sys, dt, nsteps::Int64, engine::Runners.Engine, calc::Runners.Calculation)
    pos = MolecularUtilities.vecofvec2matrix(position(sys))
    vel = MolecularUtilities.vecofvec2matrix(velocity(sys))
    inv_mass = diagm(1.0 ./ atomic_mass(sys))
    mass = diagm(atomic_mass(sys))

    phase_space = Vector{typeof(sys)}(undef, nsteps+1)
    phase_space[1] = deepcopy(sys)

    if isdir(calc.tmpfolder)
        @info "Deleting old temporary folder"
        rm(calc.tmpfolder, recursive=true)
    end
    @info "Calculating acceleration for step number 0"
    (natoms, pe, force), time_taken, _, _, _ = @timed Runners.execute(engine, calc, sys)
    energies = Matrix{typeof(pe)}(undef, nsteps+1, 3)
    energies[1, 1] = 0.5 * tr(vel * mass * vel')
    energies[1, 2] = pe
    energies[1, 3] = energies[1,1] + energies[1,2]
    @info "Step number 0 took $(round(time_taken; digits=3)) sec. KE = $(energies[1,1]). PE = $(energies[1,2]). Total Energy = $(energies[1,3])."
    acceleration = force * inv_mass
    for j = 1:nsteps
        vel .+= acceleration * dt / 2
        pos .+= vel * dt
        position(sys) .= MolecularUtilities.matrix2vecofvec(pos)
        @info "Calculating acceleration for step number $j"
        (natoms, pe, force), time_taken, _, _, _ = @timed Runners.execute(engine, calc, sys)
        acceleration = force * inv_mass
        vel .+= acceleration * dt / 2
        energies[j+1, 1] = 0.5 * tr(vel * mass * vel')
        energies[j+1, 2] = pe
        energies[j+1, 3] = energies[j+1,1] + energies[j+1,2]
        @info "Step number $j took $(round(time_taken; digits=3)) sec. KE = $(energies[j+1,1]). PE = $(energies[j+1,2]). Total Energy = $(energies[j+1,3])."
        velocity(sys) .= MolecularUtilities.matrix2vecofvec(vel)
        phase_space[j+1] = deepcopy(sys)
    end
    phase_space, energies
end

end