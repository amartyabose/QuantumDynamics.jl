using PackageCompiler

default_compile_dir() = joinpath(homedir(), ".julia", "sysimages")

default_compile_filename() = "sys_qd.so"

default_compile_path() = joinpath(default_compile_dir(), default_compile_filename())

function compile_note(; dir=default_compile_dir(), filename=default_compile_filename())
    path = joinpath(dir, filename)
    return """
    You will be able to start Julia with a compiled version of QuantumDynamics using:

    ```
    ~ julia --sysimage $path
    ```

    and you should see that the startup times and JIT compilation times are substantially improved when you are using QuantumDynamics.

    In unix, you can create an alias with the Bash command:

    ```
    ~ alias julia_qd="julia --sysimage $path -e 'using QuantumDynamics' -i"
    ```

    which you can put in your `~/.bashrc`, `~/.zshrc`, etc. This also executes
    `using QuantumDynamics` so that QuantumDynamics is loaded and ready to use, you can leave off `
    -e 'using QuantumDynamics' -i` if you don't want that. Then you can start Julia with a
    version of QuantumDynamics installed with the command:

    ```
    ~ julia_qd
    ```

    Note that if you update QuantumDynamics to a new version, for example with `using
    Pkg; Pkg.update("QuantumDynamics")`, you will need to run the `QuantumDynamics.compile()`
    command again to recompile the new version of QuantumDynamics.
    """
end

function compile(;
    dir::AbstractString=default_compile_dir(),
    filename::AbstractString=default_compile_filename()
)
    if !isdir(dir)
        println("""The directory "$dir" doesn't exist yet, creating it now.""")
        println()
        mkdir(dir)
    end
    path = joinpath(dir, filename)
    println(
        """Creating the system image "$path" containing the compiled version of QuantumDynamics. This may take a few minutes.""",
    )
    @show @__DIR__
    create_sysimage(
        :QuantumDynamics;
        sysimage_path=path,
        precompile_execution_file=joinpath(@__DIR__, "snoop/snoop.jl")
    )
    println(compile_note(; dir=dir, filename=filename))
    return path
end

@doc """
    QuantumDynamics.compile(; dir = "$(default_compile_dir())",
                       filename = "$(default_compile_filename())")

Compile QuantumDynamics.jl with [PackageCompiler](https://julialang.github.io/PackageCompiler.jl/dev/).
This will take some time, perhaps a few minutes.

This will create a system image containing the compiled version of QuantumDynamics
located at `dir/filename`, by default `$(default_compile_path())`.

$(compile_note())
""" compile