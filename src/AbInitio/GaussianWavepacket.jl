module GaussianWavepacket

using LinearAlgebra
using AtomsIO
using Unitful, UnitfulAtomic

using ..Runners, ..MolecularUtilities

abstract type HessianCalculator end
struct SingleHessian <: HessianCalculator
    hessian::AbstractMatrix
end
evaluate(ch::SingleHessian, _::AtomsIO.ExtXYZ.Atoms) = ch.hessian

function propagate_Ps_Qs(P, Q, traj::AbstractVector{<:AtomsIO.ExtXYZ.Atoms}, hess::HessianCalculator, numsteps::Int64, dt)
    Ps = Array{typeof(P[1,1])}(undef, size(P, 1), size(P, 2), numsteps+1)
    Qs = Array{typeof(Q[1,1])}(undef, size(Q, 1), size(Q, 2), numsteps+1)
    Ps[:, :, 1] = copy(P)
    Qs[:, :, 1] = copy(Q)
    As = Array{typeof((P * inv(Q))[1,1])}(undef, size(P, 1), size(Q, 1), numsteps+1)
    As[:, :, 1] = -1im * Ps[:,:,1] * inv(Qs[:,:,1])

    for j = 1:numsteps
        acc = - evaluate(hess, traj[j]) * Qs[:,:,j]
        Ps[:, :, j+1] = Ps[:, :, j] + acc * dt/2
        Qs[:, :, j+1] = Qs[:, :, j] + Ps[:,:,j+1] * dt
        acc = - evaluate(hess, traj[j]) * Qs[:, :, j+1]
        Ps[:, :, j+1] += acc * dt/2
        As[:, :, j+1] = -1im * Ps[:,:,j+1] * inv(Qs[:,:,j+1])
    end
    As, Ps, Qs
end

function get_momentum(traj::AbstractVector{<:AtomsIO.ExtXYZ.Atoms}, vecs)
    sqrtmass = diagm(sqrt.(atomic_mass(traj[1])))
    p = Matrix{typeof(velocity(traj[1])[1][1] * sqrtmass[1,1])}(undef, size(vecs, 2), length(traj))
    for (j, sys) in enumerate(traj)
        vel = MolecularUtilities.vecofvec2matrix(velocity(sys))
        ormat = MolecularUtilities.get_rotation_matrix(; sys, ref=traj[1])
        val = vecs' * reshape(ormat * vel * sqrtmass, length(vel))
        p[:, j] .= val
    end
    p
end

function get_corr(traj::AbstractVector{<:AtomsIO.ExtXYZ.Atoms}, dt, hess::HessianCalculator, A0, init_eigvec)
    p = get_momentum(traj, init_eigvec)
    Q0 = inv(sqrt.(imag.(A0))) * (1.0 + 0.0im)
    P0 = A0 * Q0
    nsteps = length(traj) - 1
    npoints = length(traj)
    As, Ps, Qs = propagate_Ps_Qs(P0, Q0, traj, hess, nsteps, dt)
    lagrangian = [0.5*p[:, j]'*p[:, j] - (traj[j].system_data.PE * 1u"Eh_au" - traj[1].system_data.PE * 1u"Eh_au") for j=1:npoints]
    S = 0u"Eh_au*fs"
    corr = zeros(ComplexF64, npoints)
    for j in axes(As, 3)
        v, _ = eigen(ustrip.(As[:,:,j]))
        corr[j] =  exp(-1im * p[:, j]' * Qs[:, :, j] * inv(Ps[:, :, j]) * p[:, j] / Unitful.ħ + 2im * S / Unitful.ħ) / (prod(sqrt.(v)) * det(ustrip.(Qs[:,:,j])))
        if j<npoints
            S += (lagrangian[j] + lagrangian[j+1]) * dt / 2
        end
    end
    0.0u"fs":2dt:2dt*nsteps, corr
end

function get_corr_single_hessian(engine::Runners.Engine, sys::AtomsIO.ExtXYZ.Atoms, traj::AbstractVector{<:AtomsIO.ExtXYZ.Atoms}, dt, init_hess_file::String)
    mwhess, vals, vecs = MolecularUtilities.get_hessian_normal_modes(engine, sys, init_hess_file)
    A0 = diagm(1im * sqrt.(vals))
    hess = SingleHessian(vecs' * mwhess * vecs)
    ts, corr = get_corr(traj, dt, hess, A0, vecs)
    denergy = (sys.system_data.PE - traj[1].system_data.PE) * 1u"Eh_au"
    ts, corr .* exp.(1im * denergy .* ts ./ Unitful.ħ)
end

end