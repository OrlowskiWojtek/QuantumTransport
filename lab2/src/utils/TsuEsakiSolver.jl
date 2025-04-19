include("TransferMatrix.jl")

struct TsuEsaki
    mu_s::Float64
    mu_d::Float64
    temp::Float64
    bias::Float64
    material::Function
end

function solve(ts::TsuEsaki)
    constants = getConstants()
    energies = collect(LinRange(1e-6, ts.mu_s, 1000))
    
    integral = 0.
    for (idx, ene) in enumerate(energies)
        qs = initQuantumSystem(ene, ts.material, applied_bias = ts.bias)
        r = solve(qs)
        nomin   = 1. + exp((ts.mu_s - ene) / (constants.BoltzmannConstant * ts.temp))
        denomin = 1. + exp((ts.mu_d - qs.applied_bias - ene) / (constants.BoltzmannConstant * ts.temp))
        integral += r.transmittance * log(nomin / denomin)
    end

    multiplier = constants.GaAsEm * constants.BoltzmannConstant * ts.temp / (2 * π^2)

    return integral * multiplier
end
