"""
    create_and_select_group(base, new_group)

Checks if `new_group` exists in `base`. Selects it if it does, else creates it.
"""
function create_and_select_group(base, new_group)
    if !haskey(base, new_group)
        create_group(base, new_group)
    end
    return base[new_group]
end

"""
    check_or_insert_value(base, variable, value)

Inserts `value` into `base[variable]`.  If it already exists, checks if the value is correct, and throws if not.
"""
function check_or_insert_value(base, variable, value)
    if !haskey(base, variable)
        base[variable] = value
    else
        @assert read_dataset(base, variable) == value "$(variable)'s value is not matching."
    end
end

function merge_HDF5(source, destination)
    for ks in keys(source)
        sks = source[ks]
        if typeof(sks) == HDF5.Dataset
            @info "HDF5.Dataset, $(ks) found. Merging."
            check_or_insert_value(destination, ks, read_dataset(source, ks))
        elseif typeof(sks) == HDF5.Group
            @info "HDF5.Group, $(ks) found. Merging."
            dgroup = create_and_select_group(destination, ks)
            merge_HDF5(sks, dgroup)
        end
    end
end

"""
    merge_into(source::String, destination::String)

Merge data from the HDF5 file at `source` to the one at `destination`.
"""
function merge_into(source::String, destination::String)
    fdestination = isfile(destination) ? h5open(destination, "r+") : h5open(destination, "w")
    fsource = h5open(source, "r")
    merge_HDF5(fsource, fdestination)
    close(fdestination)
    close(fsource)
end

function propagate_density_matrices(; filename::AbstractString, path::AbstractString, prop_name::AbstractString, init_states::Dict{<:AbstractString,<:AbstractMatrix{<:Complex}}, time_name::AbstractString="time")
    fsource = h5open(filename, "r+")
    dat_group = fsource[path]
    propagators = read_dataset(dat_group, prop_name)
    time = read_dataset(dat_group, time_name)
    ntimes = size(propagators, 1)
    ρs = []
    for (outname, ρ0) in init_states
        @info "Propagating initial state named $(outname)"
        display(ρ0)
        _, ρ = apply_propagator(; propagators, ρ0, ntimes, dt=1.0)
        if haskey(dat_group, outname)
            delete_object(dat_group, outname)
        end
        dat_group[outname] = ρ
        push!(ρs, ρ)
    end
    close(fsource)
    time, Dict(keys(init_states) .=> ρs)
end
