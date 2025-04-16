using CairoMakie
using DelimitedFiles
using ColorSchemes: colorschemes

## loading data

hxs = collect(20:20:500)

energies = Matrix{Float64}(undef, 10, 0)
for (idx, hw) in enumerate(hxs)
    filename = "../data/eigenvalues_task_4_hwx_$hw"
    row = readdlm(filename)
    energies = hcat(energies, row)
end

##
with_theme(
    theme_latexfonts()
) do

        fig = Figure();
        ax  = Axis(fig[1,1], xlabel = "ħωₓ [meV]", ylabel = "Energia stanu [meV]") 
        cm = cgrad(:managua10, 10)
        for (idx,ene) in enumerate(eachrow(energies))
            lines!(ax, hxs, ene, label = "Stan $idx", linewidth = 3, color = cm[idx])
        end

        Legend(fig[1,2], ax)

        xlims!(ax, (20, 500))
        save("plots/energies_10_states.pdf", fig)
        display(fig)
    end
