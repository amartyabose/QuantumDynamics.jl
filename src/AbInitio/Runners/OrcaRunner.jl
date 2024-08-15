using JSON3
using AtomsIO
using Unitful, UnitfulAtomic

function execute(engine::Orca, input::String, tmpfolder::String)
    original_dir = pwd()
    if !isdir(tmpfolder)
        mkdir(tmpfolder)
    end
    cd(tmpfolder)
    rm.(filter(x -> splitext(x)[2] != ".gbw", readdir(".")))
    open("input.inp", "w") do io
        write(io, input)
    end
    run(pipeline(`$(engine.executable) input.inp`, stdout="input.log", stderr="input.err"))
    # run(pipeline(`orca_2json input -property`, stdout="orca2json.log", stderr="orca2json.err"))
    props = filter(x->!startswith(x, "#"), split(read("input.engrad", String), "\n"))
    cd(original_dir)
    natoms = parse(Int64, props[1])
    energy = parse(Float64, props[2]) * 1u"Eh_au"
    force = -reshape(parse.(Float64, props[3:2+3*natoms]), 3, natoms) .* 1u"Eh_au/a0_au"
    natoms, energy, force
end
function execute(engine::Orca, calc::Calculation, system::AtomsIO.ExtXYZ.Atoms)
    input = get_input(engine, calc, system)
    execute(engine, input, calc.tmpfolder)
end

function get_input(::Orca, calc::DFT, system::AtomsIO.ExtXYZ.Atoms)
    input = "! $(calc.method) $(calc.basis)\n"
    input *= calc.numerical_gradient ? "! engrad numgrad\n\n" : "! engrad\n\n"
    if calc.num_procs > 1
        input *= "%pal nprocs $(calc.num_procs) end\n\n"
    end
    if calc.excited_state
        input *= """
%$(calc.excited_state_method)
    nroots $(calc.num_roots)
    iroot  $(calc.iroot)
end

"""
    end
    input *= "* xyz 0 1\n"
    atom_names = String.(system.atom_data.atomic_symbol)
    positions = ustrip.(system.atom_data.position)
    for (atom, pos) in zip(atom_names, positions)
        input *= "$atom"
        for p in pos
            input *= " $p"
        end
        input *= "\n"
    end
    input * "*"
end