module Utilities

function trapezoid(x, y; discrete::Bool=false)
    if discrete
        return sum(y)
    end
    sumvar = zero(y[1])
    for (a, b) in zip(y[2:end], y)
        sumvar += a + b
    end
    sumvar / 2 * (x[2] - x[1])
end

"""
    unhash_path(path_num::Int, ntimes::Int, sdim::Int)
Construct a path for a system with `sdim` dimensions, corresponding to the number `path_num`, with `ntimes` time steps.
"""
function unhash_path(path_num::Int, ntimes::Int, sdim::Int)
    path_num -= 1
    states = zeros(UInt8, ntimes + 1)
    for j in 1:ntimes+1
        @inbounds states[j] = path_num % sdim
        path_num = path_num ÷ sdim
    end
    states .+ 1
end

"""
ExternalField provides an abstract interface for encoding an external field, `V(t)`, interacting with the system through the operator, `coupling_op`.
"""
struct ExternalField
    V :: Function
    coupling_op :: Matrix{ComplexF64}
end

end