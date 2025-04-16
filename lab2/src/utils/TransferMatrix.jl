function TransferMatrix(qs::QuantumSystem,
                        n::Int64)
    T::Matrix{Complex{Float64}} = Matrix{Complex{Float64}}(undef, 2, 2)
    T[1, 1] = 0.5 * (1 + (qs.k[n + 1] * qs.m[n]) / (qs.k[n] * qs.m[n + 1])) * exp(im * (qs.k[n + 1] - qs.k[n]) * qs.x[n])
    T[1, 2] = 0.5 * (1 - (qs.k[n + 1] * qs.m[n]) / (qs.k[n] * qs.m[n + 1])) * exp(-im * (qs.k[n + 1] + qs.k[n]) * qs.x[n]) 
    T[2, 1] = 0.5 * (1 - (qs.k[n + 1] * qs.m[n]) / (qs.k[n] * qs.m[n + 1])) * exp(im * (qs.k[n + 1] + qs.k[n]) * qs.x[n]) 
    T[2, 2] = 0.5 * (1 + (qs.k[n + 1] * qs.m[n]) / (qs.k[n] * qs.m[n + 1])) * exp(-im * (qs.k[n + 1] - qs.k[n]) * qs.x[n])
    return T
end

include("SystemResults.jl")
include("QuantumSystem.jl")

function solve(qs::QuantumSystem)
    results = QSResults(0.,0.)
    matrix = TransferMatrix(qs, 1)

    for i in 2:length(qs.x)-1
        matrix *= TransferMatrix(qs, i)
    end
    
    trans_first = (real(qs.k[end])*qs.m[begin]) / (real(qs.k[begin]) * qs.m[end])
    trans_second = 1. / abs(matrix[1,1])^2

    results.transmittance =  trans_first * trans_second
    results.reflectance  = abs(matrix[2,1])^2 / abs(matrix[1,1])^2

    return results
end

