module Utilities

function trapezoid(x, y; discrete::Bool=false)
    if discrete
        return sum(y)
    end
    sum = zero(y[1])
    for (a, b) in zip(y[2:end], y)
        sum += a + b
    end
    sum / 2 * (x[2] - x[1])
end

function unhash_path(path_num::Int, ntimes, sdim)
    path_num -= 1
    states = zeros(UInt8, ntimes + 1)
    for j in 1:ntimes+1
        @inbounds states[j] = path_num % sdim
        path_num = path_num ÷ sdim
    end
    states .+ 1
end

end