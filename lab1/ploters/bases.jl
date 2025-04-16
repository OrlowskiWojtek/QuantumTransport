using DelimitedFiles, CairoMakie

filename = "../data/basis_k_0.dat"
filename1 = "../data/basis_k_8.dat"
filename2 = "../data/basis_k_9.dat"

data = readdlm(filename, comments = true)
data1 = readdlm(filename1, comments = true)
data2 = readdlm(filename2, comments = true)

fig = Figure();
ax = Axis(fig[1,1], title = "k = 0", width = 150, height = 150, ylabel = "y [nm]");
ax1 = Axis(fig[1,2], title = "k = 8", width = 150, height = 150);
ax2 = Axis(fig[1,3], title = "k = 9", width = 150, height = 150);

Label(fig[2,1:3], "x [nm]")

x = collect(LinRange(-75.6144*0.059, -75.6144*0.059, 101))

cm = heatmap!(ax, data, interpolate = false);
heatmap!(ax1, data1, interpolate = false);
heatmap!(ax2, data2, interpolate = false);

Colorbar(fig[3,1:3], cm, vertical = false)

resize_to_layout!(fig)
save("plots/basis_funcs.pdf", fig)
display(fig)

