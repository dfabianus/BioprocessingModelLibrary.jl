using Revise
using ModelingToolkit
using OrdinaryDiffEq
using BioprocessingModelLibrary.Refolding

@named kiefhaber_model = Kiefhaber()
sys = structural_simplify(kiefhaber_model)
kiefhaber_prob = ODEProblem(sys, 
[kiefhaber_model.I=>1.0, kiefhaber_model.N=>0.0, kiefhaber_model.A=>0.0], (0.0, 10.0))