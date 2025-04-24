include("utils/QuantumSystem.jl")
include("utils/Plotters.jl")
using ProgressMeter

struct TaskManager
    task1::Function
    task2::Function
    task3::Function
end

function task1( disp = false)
    constants = getConstants()
    energies = collect(LinRange(0.,1.,100))

    # task 1.1
    transmittance = Float64[]
    reflectance = Float64[]

    function material(pos)
        if((pos > 5. / constants.BohrRadius) && (pos < 10. / constants.BohrRadius))
            return GaAlAs::Material
        else
            return GaAs::Material
        end
    end

    for ene in energies
        ene = ene / 27.2116
        qs = initQuantumSystem(ene, material, change_mass = false)
        r = solve(qs)
        push!(transmittance, r.transmittance)
        push!(reflectance, r.reflectance)
    end

    f = plotTransRefl(energies, transmittance, reflectance, disp);
    save("../plots/single_barrier_const_mass.pdf", f)

    # task 1.2
    transmittance = Float64[]
    reflectance = Float64[]

    for ene in energies
        ene = ene / 27.2116
        qs = initQuantumSystem(ene, material)
        r = solve(qs)
        push!(transmittance, r.transmittance)
        push!(reflectance, r.reflectance)
    end

    f = plotTransRefl(energies, transmittance, reflectance, disp);
    save("../plots/single_barrier.pdf", f)
end

include("utils/TsuEsakiSolver.jl")

function task2(disp)
    constants = getConstants()
    energies = collect(LinRange(0.,1.,1000))

    # task 2.1
    transmittance = Float64[]
    reflectance = Float64[]

    function material(pos)
        if((pos > 1. / constants.BohrRadius) && (pos < 6. / constants.BohrRadius))
            return GaAlAs::Material
        elseif ((pos > 9. / constants.BohrRadius) && pos < 14. / constants.BohrRadius)
            return GaAlAs::Material
        else
            return GaAs::Material
        end
    end

    for ene in energies
        ene = ene / 27.2116
        qs = initQuantumSystem(ene, material)
        r = solve(qs)
        push!(transmittance, r.transmittance)
        push!(reflectance, r.reflectance)
    end

    f = plotTransRefl(energies, transmittance, reflectance, disp);
    save("../plots/double_barrier.pdf", f)

    # task 2.2
    transmittance = Float64[]
    reflectance = Float64[]

    function material(pos)
        if((pos > 1. / constants.BohrRadius) && (pos < 6. / constants.BohrRadius))
            return GaAlAs::Material
        elseif ((pos > 9. / constants.BohrRadius) && pos < 14. / constants.BohrRadius)
            return GaAlAs::Material
        else
            return GaAs::Material
        end
    end

    qs = initQuantumSystem(0.0, material, applied_bias = 50. / constants.HartreeTomV)
    f = plotQuantumSystem(qs, disp) # checking applied bias function
    save("../plots/double_barrier_with_bias.pdf", f)

    # for defined system solve TsuEsakiFormula
    V_bias = collect(LinRange(0, 500, 200))
    currents = Vector{Float64}(undef, length(V_bias))
    mu_s = 87. / constants.HartreeTomV
    mu_d = 87. / constants.HartreeTomV
    temp = 77.

    @showprogress desc = "Computing from Tsu-Esaki formula" for (idx, V) in enumerate(V_bias)
        ts = TsuEsaki(mu_s, mu_d, temp, V / constants.HartreeTomV, material)
        currents[idx] = solve(ts)
    end

    f = plotCurrentVoltage(currents, V_bias, disp)
    save("../plots/iv_characteristic.pdf", f)
end

include("utils/QuantumPointContact.jl")

function task3(disp)
    qpc = initQPC()

    f = plotQPCPotential(qpc, disp)
    save("../plots/qpc_potential_map.pdf", f)

    energies = calculateEffectivePotential(5, qpc)
    f = plotQPCEnergies(energies, qpc, disp)
    save("../plots/qpc_states_energies.pdf", f)

    ene, cond = calcConductance(qpc)
    f = plotConductances(ene, cond, true)
    save("../plots/qpc_conductance.pdf", f)

    # task 3.3
    c = getConstants()
    E1 = 50   / c.HartreeTomV
    E2 = 100  / c.HartreeTomV

    voltages = LinRange(0, 25, 100)
    cond_e1 = Vector{Float64}(undef, length(voltages))
    cond_e2 = Vector{Float64}(undef, length(voltages))
    nstates = 5

    @showprogress desc = "Calculating conductance vs Vg" for (idx, vgate) in enumerate(voltages)
        qpc = initQPC(vgate * 1000)
        ce1 = getConductance(E1, nstates, qpc, 100)
        ce2 = getConductance(E2, nstates, qpc, 100)
        cond_e1[idx] = ce1
        cond_e2[idx] = ce2
    end
    
    f = plotConductanceVsVoltage(voltages, cond_e1, cond_e2, disp)
    save("../plots/cond_vs_voltage.pdf", f)
end

function getTaskManager()

    return TaskManager(task1, task2, task3)
end
