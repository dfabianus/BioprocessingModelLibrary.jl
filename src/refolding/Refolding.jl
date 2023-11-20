module Refolding

export Kiefhaber, Kiefhaber_network

using ModelingToolkit
using Catalyst

@parameters t
Dt = Differential(t)

include("kiefhaber_1991.jl")

end
