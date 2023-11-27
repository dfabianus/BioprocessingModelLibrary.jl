module Refolding

export Kiefhaber, Kiefhaber_B, Kiefhaber_network
export INAM, INAM_network
export FLUMO_COMBINED_3, FLUMO_COMBINED_4

using ModelingToolkit
using Catalyst

@parameters t
Dt = Differential(t)

include("kiefhaber_1991.jl")
include("INA+M.jl")
include("flumo_models.jl")

end
