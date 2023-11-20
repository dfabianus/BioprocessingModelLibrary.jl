using Revise
using ModelingToolkit
using OrdinaryDiffEq
using BioprocessingModelLibrary.Refolding
using Plots

@named kiefhaber_model = Kiefhaber()
sys = structural_simplify(kiefhaber_model)
kiefhaber_prob = ODEProblem(sys, 
[sys.I=>1.0, sys.N=>0.0, sys.A=>0.0], (0.0, 10.0))
sol = solve(kiefhaber_prob, Tsit5(), saveat=0.1)
