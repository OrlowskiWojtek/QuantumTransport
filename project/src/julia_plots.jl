using CairoMakie
using DelimitedFiles

set_theme!(merge(theme_dark(), theme_latexfonts()))
##
# do przepisania ręcznie:

Emax    = 0.7
Emin    = -0.7
En      = 30
enes = LinRange(Emin, Emax, En)
V_n = 100
Vmin = -0.1
Vmax = 0.1
V_np_vec = LinRange(Vmin, Vmax, V_n)
selected_line_idx = 2

cond = readdlm("data/conductance_v2.dat")

fig = Figure();
ax = Axis(fig[1,1], title = "G", xlabel = L"V_{np} [\text{eV}]", ylabel = L"E [\text{eV}]");
ax_diff = Axis(fig[1,2], title = L"\frac{\partial G}{\partial V_{np}}", xlabel = L"V_{np} [\text{eV}]")
ax_diff_lines = Axis(fig[3,1:2], title = L"\frac{\partial G}{\partial V_{np}}", xlabel = L"V_{np} [\text{eV}]")

dV = V_np_vec[2] - V_np_vec[1]
cm = heatmap!(ax, V_np_vec, enes, cond, colormap = :bluesreds);
cm_diff = heatmap!(ax_diff, V_np_vec[begin:end-1], enes, diff(cond, dims = 1) ./ dV, colormap = :bluesreds);
hlines!(ax_diff, [enes[begin+selected_line_idx]],linestyle = :dash, color = :black)

selected_line = diff(cond, dims = 1)[:, selected_line_idx]
lines!(ax_diff_lines, V_np_vec[begin:end-1], selected_line)


Colorbar(fig[2, 1], cm, vertical = false)
Colorbar(fig[2, 2], cm_diff, vertical = false)

display(fig)
save("../figures/conductance_long.pdf", fig)

## changing magnetic field plot

Ene = 1;
cond = readdlm("data/cond_b_temp.dat")
V_n = 100
V_min = -1.
V_max = 1.
V_np_vec = collect(LinRange(V_min, V_max, V_n))

B_n = 100
B_min = -1.
B_max = 1.
B_vec = LinRange(B_min, B_max, B_n)

fig = Figure();
ax = Axis(fig[1,1], aspect = 1., title = "G", xlabel = L"V_{np} [\text{eV}]", ylabel = L"B [\text{T}]");

cm = heatmap!(ax, V_np_vec, B_vec, cond, colormap = :bluesreds);

Colorbar(fig[2, 1], cm, vertical = false)

#save("../figures/conductance_magnetic.pdf", fig)
display(fig);
