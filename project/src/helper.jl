using CairoMakie

##
L = 200 / 0.0529
W = 100 / 0.0529
d = 5. / 0.0529

xs = LinRange(-L/2, L/2, 200);

pot(x) =@. 5 * (tanh((x + L/4)/d) - tanh((x - L/4)/d)) 

p = pot(xs)

lines(xs, p)
