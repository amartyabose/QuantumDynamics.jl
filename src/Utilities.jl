module Utilities

function unhash_path(path_num::Int, ntimes, sdim)
    path_num -= 1
    states = zeros(UInt8, ntimes+1)
    for j in 1:ntimes+1
        states[j] = path_num % sdim
        path_num = path_num ÷ sdim
    end
    states .+ 1
end

end