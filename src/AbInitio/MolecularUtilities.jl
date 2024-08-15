module MolecularUtilities

using ..Runners

using LinearAlgebra
using AtomsIO
using Unitful, UnitfulAtomic

function vecofvec2matrix(x)
    ans = Matrix{eltype(x[1])}(undef, length(x[1]), length(x))
    for (j, v) in enumerate(x)
        ans[:, j] .= v
    end
    ans
end

matrix2vecofvec(x) = [x[:, j] for j in axes(x, 2)]

get_CoM(mass, pos) = sum([m*r for (m, r) in zip(mass, pos)]) / sum(mass)

function get_CoM(sys::AtomsIO.ExtXYZ.Atoms)
    mass = atomic_mass(sys)
    pos = position(sys)
    get_CoM(mass, pos)
end

function shift_system!(sys::AtomsIO.ExtXYZ.Atoms, new_origin)
    newpos = [pos .- new_origin for pos in position(deepcopy(sys))]
    position(sys) .= newpos
    nothing
end

function get_rotation_matrix(; sys::AtomsIO.ExtXYZ.Atoms, ref::AtomsIO.ExtXYZ.Atoms)
    shift_system!(sys, get_CoM(sys))
    shift_system!(ref, get_CoM(ref))
    pos_sys = vecofvec2matrix(position(sys))
    pos_ref = vecofvec2matrix(position(ref))
    covar = ustrip.(pos_sys * transpose(pos_ref))
    U, _, V = svd(covar)
    @info det(U), det(V)
    V * transpose(U)
end

function get_hessian_normal_modes(::Runners.Orca, sys::AtomsIO.ExtXYZ.Atoms, hessian_file::String)
    mass = atomic_mass(sys)
    natoms = length(mass)
    hess_unit = 1u"Eh_au" / 1u"a0_au"^2
    hessian_size = 3*natoms
    hess = zeros(typeof(hess_unit), hessian_size, hessian_size)

    hess_content = filter(x->length(x)!=0, split(read(hessian_file, String), "\n"))
    hess_line = findfirst(x->x=="\$hessian", hess_content)

    ncols = 5

    for (sn, k) in enumerate(1:5:hessian_size)
        for j = 1:hessian_size
            hess[j, k:min(ncols-1+k, hessian_size)] = parse.(Float64, split(hess_content[hess_line+2+j+(sn-1)*(hessian_size+1)])[2:end]) .* 1u"Eh_au*a0_au^-2"
        end
    end

    mass_weighted_hess = zeros(typeof(hess_unit / mass[1]), hessian_size, hessian_size)
    for j=1:hessian_size, k=1:hessian_size
        mass_weighted_hess[j,k] = hess[j,k] / sqrt(mass[(j-1)_3+1] * mass[(k-1)_3 + 1])
    end

    vib_freq_row_num = findfirst(x -> x=="\$vibrational_frequencies", hess_content)

    freq_vals = filter(x->x==0.0u"eV", parse.(Float64, [split(hs)[2] for hs in hess_content[vib_freq_row_num+2:vib_freq_row_num+hessian_size+1]]) .* 1u"cm^-1" .* Unitful.h .* Unitful.c)
    vals, vecs = eigen(ustrip.(mass_weighted_hess))
    vals = vals[length(freq_vals)+1:end] .* unit(mass_weighted_hess[1,1])
    vecs = vecs[:, length(freq_vals)+1:end]

    mass_weighted_hess, vals, vecs
end

end