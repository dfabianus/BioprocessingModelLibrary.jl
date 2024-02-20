module Upstream

export ColiModel

using ModelingToolkit
using Catalyst

@parameters t
Dt = Differential(t)

include("colimodel_A.jl")

end