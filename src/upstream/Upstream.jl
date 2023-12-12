module Upstream

export Pablymueller

using ModelingToolkit
using Catalyst

@parameters t
Dt = Differential(t)

include("pablymueller.jl")

end