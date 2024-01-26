module Refolding

export Kiefhaber, Kiefhaber_B, Kiefhaber_network
export INAM, INAM_network
export FLUMO_COMBINED_3, FLUMO_COMBINED_4, FLUMO_COMBINED_5, FLUMO_COMBINED_6
export FLUMO_pulsing_events
export FLUMO_MECH
export FLUMO_SOFTSENSOR
export CIGSTEADY

using ModelingToolkit
using Catalyst
import IfElse

@parameters t
Dt = Differential(t)

include("kiefhaber_1991.jl")
include("INA+M.jl")
include("flumo_models.jl")

end