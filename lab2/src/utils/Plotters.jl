include("QuantumSystem.jl")
include("TransferMatrix.jl")

using CairoMakie

set_theme!(theme_latexfonts())

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
    ax = Axis(fig[1,1], xlabel = "Vbias [meV]", ylabel = "j")
    lines!(voltage, current)

    disp && display(fig)

    return fig
end

include("QuantumPointContact.jl")

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
    ax = Axis(fig[1,1], xlabel = "x [nm]", ylabel = "Energy [meV]")
    c = getConstants()

    xvec = collect(LinRange(0, qpc.length, size(energies, 2))) .* c.BohrRadius
    for (idx, state) in enumerate(eachrow(energies))
        lines!(ax, xvec, state .* c.HartreeTomV, label = "State $idx")
    end

    axislegend()
    disp && display(fig)

    return fig
end

function plotConductances(energies, conductances, disp = false)::Figure
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "Energy [eV]", ylabel = "G [2e²/h]")
    c = getConstants()

    lines!(ax, energies .* c.HartreeTomV, conductances)
    
    disp && display(fig)

    return fig
end

function plotConductanceVsVoltage(voltage, cond_e1, cond_e2, disp = false)::Figure
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "Voltage [V]", ylabel = "G [2e²/h]")

    lines!(ax, voltage, cond_e1, label = "E = 50 meV")
    lines!(ax, voltage, cond_e2, label = "E = 100 meV")
    
    axislegend()

    disp && display(fig)
    return fig
end
