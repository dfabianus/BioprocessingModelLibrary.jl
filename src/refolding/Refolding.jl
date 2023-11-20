module Refolding

export Kiefhaber

using ModelingToolkit

@parameters t
Dt = Differential(t)

include("kiefhaber_1991.jl")

end
