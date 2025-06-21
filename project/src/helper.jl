using CairoMakie

##
L = 200 / 0.0529
W = 100 / 0.0529
d = 5. / 0.0529

xs = LinRange(-L/2, L/2, 200);

pot(x) =@. 5 * (tanh((x + L/8)/d) - tanh((x - L/8)/d)) / 2.

p = pot(xs)

fig = Figure();
ax = Axis(fig[1,1], yticks = ([5.], [L"V_{np}"]), xlabel = "x [nm]", ylabel = "potential [eV]");

lines!(ax, xs * 0.0529, p)

display(fig)
save("../figures/barrier_d5.pdf", fig)
