include("utils/TransferMatrix.jl")
##
transmittance = Float64[]
reflectance = Float64[]
energies = collect(LinRange(0.,1.,100))

for ene in energies
    ene = ene / 27.2116
    qs = initQuantumSystem(ene)
    r = solve(qs)
    push!(transmittance, r.transmittance)
    push!(reflectance, r.reflectance)
end

using CairoMakie

fig = Figure();
ax = Axis(fig[1,1])

lines!(ax, energies, transmittance, color = :red, label = "transmittance")
lines!(ax, energies, reflectance, color = :blue, label = "reflectance")
axislegend()

display(fig)



