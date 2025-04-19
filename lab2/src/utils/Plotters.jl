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
