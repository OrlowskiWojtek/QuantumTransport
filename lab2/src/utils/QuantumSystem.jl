include("constants.jl")

@enum Material GaAs GaAlAs

struct QuantumSystem
    U::Vector{Float64}
    k::Vector{ComplexF64}
    x::Vector{Float64}
    m::Vector{Float64}
    E::Float64
    constants::NamedTuple
    material::Function # function that returns material enum
    xmin::Float64
    xmax::Float64 
    applied_bias::Float64
end

function genPotential(x::Vector{Float64}, params, material; bias::Function)
    V = Vector{Float64}(undef, length(x))
    for (idx, xval) in enumerate(x)
        if (material(xval) == GaAs::Material)
            V[idx] = 0. - bias(xval)
        elseif (material(xval) == GaAlAs::Material)
            V[idx] = params.V - bias(xval)
        end
    end

    return V
end

function genMass(x::Vector{Float64}, params, material; change_mass = true)
    M = Vector{Float64}(undef, length(x))
    for (idx, xval) in enumerate(x)
        if (material(xval) == GaAs::Material)
            M[idx] = params.GaAsEm
        elseif (material(xval) == GaAlAs::Material)
            M[idx] = params.GaAlAsEm
            if(!change_mass) M[idx] = params.GaAsEm end
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

function initQuantumSystem(Energy::Float64, material::Function; change_mass = true, applied_bias::Float64 = 0.)
    constants = getConstants()
    xmin = 0. / constants.BohrRadius
    xmax = 15. / constants.BohrRadius
    N = 1000

    bias(x) = applied_bias * x / xmax
    x = collect(LinRange(xmin, xmax , N))
    V = genPotential(x, constants, material, bias = bias)
    mass = genMass(x, constants, material, change_mass = change_mass)
    k = genKvector(V, mass, Energy, constants)

    qs = QuantumSystem(V, k, x, mass, Energy, constants, material, xmin, xmax, applied_bias)
end
