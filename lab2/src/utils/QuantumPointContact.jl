include("constants.jl")

using LinearAlgebra

struct QuantumPointContact
    ε::Float64  # permitivity
    d::Float64  # distance between the gates
    Vg::Float64 # gate voltage
    width::Float64
    length::Float64 
    gate_begin::Float64 # relative to length
    gate_end::Float64   # relative to length
    gate_spacing::Float64 # relative to width
end

function f(u, v, qpc::QuantumPointContact)
    return 1. / (2π * qpc.ε) * atan( u*v / (qpc.d * sqrt(qpc.d^2 + u^2 + v^2)) )
end

function V(x, y, qpc::QuantumPointContact)
    l = qpc.length * qpc.gate_begin
    r = qpc.length * qpc.gate_end

    # gate one -> bottom one
    t = - qpc.gate_spacing * qpc.width / 2. 
    b = - qpc.width * 2 # overrate this one | dla prowadzącego: Powinno to być napisane w instrukcji <- bardzo niejasne

    f1 = f(x - l, y - b, qpc)
    f2 = f(x - l, t - y, qpc)
    f3 = f(r - x, y - b, qpc)
    f4 = f(r - x, t - y, qpc)
    Vg1 = qpc.Vg 
    
    g1 = Vg1 * (f1 + f2 + f3 + f4)

    #gate two -> upper one
    b = qpc.gate_spacing * qpc.width / 2. 
    t = qpc.width * 2 # overrate this one

    f1 = f(x - l, y - b, qpc)
    f2 = f(x - l, t - y, qpc)
    f3 = f(r - x, y - b, qpc)
    f4 = f(r - x, t - y, qpc)
    Vg2 = qpc.Vg

    g2 = Vg2 * (f1 + f2 + f3 + f4)
    
    return g1 + g2
end

function initQPC()::QuantumPointContact
    c = getConstants()

    eps = 13.6
    d = 3. / c.BohrRadius
    Vg = 4000. / c.HartreeTomV
    W = 50. / c.BohrRadius
    L = 100. / c.BohrRadius
    gate_begin   = 0.3
    gate_end     = 0.7
    gate_spacing = 0.6

    return QuantumPointContact(eps, d, Vg, W, L, gate_begin, gate_end, gate_spacing)
end

function calculateEffectivePotential(n::Int64, qpc::QuantumPointContact; N::Int64 = 100)
    xs = collect(LinRange(0, qpc.length, N))
    ys = collect(LinRange(-qpc.width/2, qpc.width/2, N))
    dy = ys[2] - ys[1]

    c = getConstants()
    α = 1 / (2 * c.GaAsEm * dy^2)
    H = zeros(N, N)
    energies = Matrix{Float64}(undef, n, N)

    for (x_idx, x) in enumerate(xs)
        for (idx, y) in enumerate(ys)
            if(idx == 1)
                H[idx, idx] = 2 * α + V(x, y, qpc)
                H[idx, idx + 1] =  -α
            elseif(idx == length(ys))
                H[idx, idx] = 2 * α + V(x, y, qpc)
                H[idx, idx - 1] =  -α
            else       
                H[idx, idx - 1] =  -α
                H[idx, idx + 1] =  -α
                H[idx, idx] = 2 * α + V(x, y, qpc)
            end
        end
        
        x_enes = eigvals(H)
        energies[:, x_idx] = x_enes[begin:n]
    end

    return energies
end

include("TransferMatrix.jl")

function calcConductance(qpc::QuantumPointContact)
    c = getConstants()
    N = 100
    energies = LinRange(0.001 / c.HartreeTomV, 200 / c.HartreeTomV, N)
    conductances = Vector{Float64}(undef, N)

    states = 5
    for (idx, ene) in enumerate(energies)
        conductances[idx] = getConductance(ene, states, qpc, 200)
    end

    return energies, conductances
end

function getConductance(energy, states, qpc::QuantumPointContact, N::Int64)
    eff_potential = calculateEffectivePotential(states, qpc, N = N)

    c = getConstants()

    x = collect(LinRange(0., qpc.length, N))
    mass = [c.GaAsEm for i in x]
    material = x -> GaAs::Material
    temp_cond = 0.

    for state in eachrow(eff_potential)
        eff_pot = Vector(state)
        k = genKvector(eff_pot, mass, energy, c)
        qs = QuantumSystem(eff_pot, k, x, mass, energy, c, material, 0., qpc.length, 0.)
        r = solve(qs)
        temp_cond += r.transmittance
    end

    return temp_cond
end
