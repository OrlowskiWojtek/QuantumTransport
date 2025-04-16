mutable struct QSResults
    transmittance::Float64
    reflectance::Float64
end

Base.show(io::IO, r::QSResults) = print("Transmitance = $(r.transmittance) \n
                                         Reflectance  = $(r.reflectance) \n")
