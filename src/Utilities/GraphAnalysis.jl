module GraphAnalysis

using LinearAlgebra

uniq_edges(H::AbstractMatrix{<:Number}) = filter(x -> x[1]<x[2], findall(x -> x≠0.0, H))

function incidence_matrix(unique_edges::Vector{CartesianIndex{2}}, N::Int64)
    M = length(unique_edges)
    B = zeros(M, N)
    for j = 1:N
        for (k, x) in enumerate(unique_edges)
            if x[1] == j
                B[k, j] = -1.0
            elseif x[2] == j
                B[k, j] = 1.0
            end
        end
    end
    B
end
incidence_matrix(unique_edges::Vector{CartesianIndex{2}}, H::AbstractMatrix{<:Number}) = incidence_matrix(unique_edges, size(H, 1))

@enum Weights identity sgnHam absHam
function get_weights(H::AbstractMatrix{<:Number}, unique_edges::Vector{CartesianIndex{2}}, wt::Weights)
    M = length(unique_edges)
    if wt == identity
        Matrix{Float64}(I, M, M)
    end
end

graph_laplacian(incidence_matrix::Matrix{Float64}; weight_matrix::AbstractMatrix{<:Number}) = - transpose(incidence_matrix) * weight_matrix * incidence_matrix

end
