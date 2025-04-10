using DelimitedFiles, CairoMakie

hfile = "../data/H_matrix"
sfile = "../data/S_matrix"

H = readdlm(hfile)
S = readdlm(sfile)

fig = Figure();
ax_h = Axis(fig[1,1]);
ax_s = Axis(fig[1,2]);
heatmap!(ax_h, H)
heatmap!(ax_s, S)

display(fig)
