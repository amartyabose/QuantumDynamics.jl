module Blip

struct State
    Δs :: Float64
    sbar :: Vector{Float64}
    forward_ind :: Vector{UInt8}
    backward_ind :: Vector{UInt8}
end

function setup_simulation(H, svec)
    states = Dict{Float64, State}()
    for sf in svec
        for sb in svec
            Δs = sf - sb
            sbar = (sf + sb) / 2.0
        end
    end
end

end