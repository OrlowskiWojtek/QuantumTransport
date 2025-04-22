include("QuantumSystem.jl")
include("TransferMatrix.jl")

using CairoMakie

function plotQuantumSystem(qd::QuantumSystem, disp = false)::Figure
    constants = getConstants()

    fig = Figure();
    ax = Axis(fig[1,1], xlabel = "x [nm]", ylabel = "V [meV]");

    x = qd.x .* constants.BohrRadius
    V = qd.U .* constants.HartreeTomV

    lines!(ax, x, V)

    disp && display(fig)

    return fig
end

function plotTransRefl(energies, transmittance, reflectance, disp = false)::Figure
    fig = Figure();
    ax = Axis(fig[1,1], xlabel = "Energy [eV]")

    lines!(ax, energies, transmittance, color = :red, label = "transmittance")
    lines!(ax, energies, reflectance, color = :blue, label = "reflectance")
    axislegend()

    disp && display(fig)

    return fig
end

function plotCurrentVoltage(current, voltage, disp = false)::Figure
    fig = Figure();
    ax = Axis(fig[1,1], xlabel = "Vbias", ylabel = "j")
    lines!(voltage, current)

    disp && display(fig)

    return fig
end

include("QuantumPointContact.jl")

function plotQPC(qpc::QuantumPointContact, disp = false)::Figure
    fig = Figure();
    ax = Axis(fig[1,1], xlabel = "x [nm]", ylabel = "y [nm]")

    return fig
end

function plotQPCPotential(qpc::QuantumPointContact, disp = false)::Figure
    N = 100
    
    c = getConstants()
    xs = LinRange(0, qpc.length, N)
    ys = LinRange(-qpc.width / 2., qpc.width / 2., N)

    Vmap = [V(x, y, qpc) for x in xs, y in ys] * c.HartreeTomV / 1000.
    
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "x [nm]", ylabel = "y [nm]")

    cm = heatmap!(ax, xs .* c.BohrRadius, ys .* c.BohrRadius, Vmap, colormap = :thermal)
    Colorbar(fig[1,2], cm, label = "V [eV]")

    contour!(ax, xs .* c.BohrRadius, ys .* c.BohrRadius, Vmap, color = :white, levels = 30)

    disp && display(fig)

    return fig
end

function plotQPCEnergies(energies, qpc::QuantumPointContact, disp = false)::Figure
    fig = Figure()
    ax = Axis(fig[1,1])
    
    for state in eachrow(energies)
        lines!(ax, state)
    end

    disp && display(fig)

    return fig
end

function plotConductances(energies, conductances, disp = false)::Figure
    fig = Figure()
    ax = Axis(fig[1,1])

    lines!(ax, energies, conductances)
    
    disp && display(fig)

    return fig
end
