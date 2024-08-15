module Runners

abstract type Engine end
abstract type Calculation end

struct Orca <: Engine
    executable :: String
end
function Orca(; exec="\$HOME/bin/orca_6_0_0/orca")
    @assert isfile(exec)
    Orca(exec)
end

struct DFT <: Calculation
    tmpfolder::String
    method::String
    basis::String
    numerical_gradient::Bool
    num_procs::Int64
    excited_state::Bool
    excited_state_method::String
    num_roots::Int64
    iroot::Int64
end

include("OrcaRunner.jl")

end