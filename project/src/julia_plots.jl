using CairoMakie
using DelimitedFiles

# do przepisania ręcznie:

Emax    = 0.5
Emin    = -0.5
En      = 100
enes = LinRange(Emin, Emax, En)
V_n = 100
Vmin = -0.2
Vmax = 0.2
V_np_vec = LinRange(Vmin, Vmax, V_n)
selected_line_idx = 30

cond = readdlm("data/conductance_long_v3_e05_v02.dat", comments = true)

fig = Figure();
ax = Axis(fig[1,1], title = L"G \left[$\frac{e^2}{h}$\right]", xlabel = L"V_{np} [\text{eV}]", ylabel = L"E [\text{eV}]");
ax_diff = Axis(fig[1,2], title = L"\frac{\partial G}{\partial V_{np}} \left[$\frac{e^2}{h \text{eV}}$\right]", xlabel = L"V_{np} [\text{eV}]")
#ax_diff_lines = Axis(fig[3,1:2], title = L"\frac{\partial G}{\partial V_{np}}", xlabel = L"V_{np} [\text{eV}]")

dV = V_np_vec[2] - V_np_vec[1]
cm = heatmap!(ax, V_np_vec, enes, cond, colormap = :bluesreds);
cm_diff = heatmap!(ax_diff, V_np_vec[begin:end-1], enes, diff(cond, dims = 1) ./ dV, colormap = :bluesreds);

f(x) = x 
lines!(ax_diff, V_np_vec[begin:end-1],f ,linestyle = :dash, color = :black)

Colorbar(fig[2, 1], cm, vertical = false)
Colorbar(fig[2, 2], cm_diff, vertical = false)

save("../figures/conductance_long_e05_v02.pdf", fig)
display(fig)

## changing magnetic field plot

cond = readdlm("data/cond_b.dat")
E = 0.01;

V_n = size(cond)[2]
V_min = 0.00
V_max = 0.06
V_np_vec = collect(LinRange(V_min, V_max, V_n))

B_n = size(cond)[1]
B_min = 0.8
B_max = 1.
B_vec = LinRange(B_min, B_max, B_n)
B_lines_idxs = [2,8]

fig = Figure();
ax = Axis(fig[1,1], title = "Transmittance", ylabel = L"V_{np} [\text{eV}]", xlabel = L"B [\text{T}]");
ax_lines = Axis(fig[1, 2], xlabel = L"V_{np} [\text{eV}]", ylabel = "Transmittance");

cm = heatmap!(ax, B_vec, V_np_vec, cond, colormap = :bluesreds);
vlines!(ax, B_vec[B_lines_idxs][1], color = :darkred, linestyle = :dash)
vlines!(ax, B_vec[B_lines_idxs][2], color = :darkgreen, linestyle = :dash)

lines!(ax_lines, V_np_vec, cond[B_lines_idxs[1], :], color = :darkred)
lines!(ax_lines, V_np_vec, cond[B_lines_idxs[2], :], color = :darkgreen)

Colorbar(fig[2, 1], cm, vertical = false)

display(fig);
save("../figures/conductance_b_field_scan_low.png", fig)
