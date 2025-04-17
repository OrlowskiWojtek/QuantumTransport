using DelimitedFiles, CairoMakie

function load_matrices_from_file(filepath::String, N::Int)
    matrices = Matrix{Float64}[]
    current_block = Float64[]

    open(filepath, "r") do io
        for line in eachline(io)
            line = strip(line)
            if isempty(line)
                if length(current_block) == N * N
                    push!(matrices, reshape(copy(current_block), N, N))
                else
                    @warn "Pominięto blok o niewłaściwym rozmiarze: $(length(current_block)) elementów"
                end
                empty!(current_block)
            else
                row_values = parse.(Float64, split(line))
                append!(current_block, row_values)
            end
        end
        if length(current_block) == N * N
            push!(matrices, reshape(copy(current_block), N, N))
        elseif !isempty(current_block)
            @warn "Ostatni blok ma zły rozmiar: $(length(current_block)) elementów"
        end
    end

    return matrices
end

wfs = load_matrices_from_file("../data/wavefunctions_task_3", 101)
energies = readdlm("../data/eigenvalues_task_3")

##

x = collect(LinRange(-75.6144*0.059, 75.6144*0.059, 101))
with_theme(theme_latexfonts()) do
    fig = Figure();
    newaxes = [Axis(fig[(i-1) % 3, Int(floor((i - 1) / 3))]) for i in 1:6]

    for (idx, wf) in enumerate(wfs)
        heatmap!(newaxes[idx], x, x, wf, colormap = :thermal)
        newaxes[idx].title = "State $(idx-1), energy: $(energies[idx])"
    end
    Label(fig[0:2,-1], "x [nm]", rotation = π/2)
    Label(fig[3,0:1], "y [nm]")

    save("plots/wavefunctions.pdf", fig)
    display(fig)
end
