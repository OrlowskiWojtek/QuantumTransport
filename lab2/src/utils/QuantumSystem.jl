include("constants.jl")

struct QuantumSystem
    U::Vector{Float64}
    k::Vector{ComplexF64}
    x::Vector{Float64}
    m::Vector{Float64}
    E::Float64
    constants::NamedTuple
end

function genPotential(x::Vector{Float64}, params)
    V = Vector{Float64}(undef, length(x))
    for (idx, xval) in enumerate(x)
        if (xval < 5. / params.BohrRadius || xval > 10 / params.BohrRadius)
            V[idx] = 0.
        else
            V[idx] = params.V 
        end
    end

    return V
end

function genMass(x::Vector{Float64}, params)
    M = Vector{Float64}(undef, length(x))
    for (idx, xval) in enumerate(x)
        if (xval < 5. / params.BohrRadius || xval > 10 / params.BohrRadius)
            M[idx] = params.GaAsEm
        else
            M[idx] = params.GaAlAsEm
        end
    end

    return M
end

function genKvector(V::Vector{Float64}, m::Vector{Float64}, E::Float64, params)
    k_vector = Vector{ComplexF64}(undef, length(V))
    for i in eachindex(V)
        v = V[i]
        if(E > v)
            k_vector[i] = sqrt(2. * m[i] * (E - v))
        else
            k_vector[i] = im * sqrt(2. * m[i] * (v - E))
        end
    end
    return k_vector
end

function initQuantumSystem(Energy::Float64)
    constants = getConstants()
    xmin = 0. / constants.BohrRadius
    xmax = 15. / constants.BohrRadius
    N = 1000

    x = collect(LinRange(xmin, xmax , N))
    V = genPotential(x, constants)
    mass = genMass(x, constants)
    k = genKvector(V, mass, Energy, constants)
    
    qs = QuantumSystem(V, k, x, mass, Energy, constants)
end
