using Revise
using ModelingToolkit
using OrdinaryDiffEq
using BioprocessingModelLibrary.Refolding
using Plots

@named kiefhaber_model = Kiefhaber()
sys = structural_simplify(kiefhaber_model)
kiefhaber_prob = ODEProblem(sys, 
[sys.I=>1.0, sys.N=>0.0, sys.A=>0.0], (0.0, 10.0),
[sys.k_n=>0.2, sys.k_a=>0.1])
sol = solve(kiefhaber_prob, Tsit5(), saveat=0.1)
plot(sol, idxs=[sys.I, sys.N, sys.A], label=["I" "N" "A"], xlabel="Time (s)", ylabel="Concentration (M)")

oprob = ODEProblem(Kiefhaber_network, 
[:I => 1.0, :N => 0.0, :A => 0.0], 
(0.0, 10.0), 
(:k_n => 0.2, :k_a => 0.1))
osol = solve(oprob, Tsit5())
plot(osol, idxs=[Kiefhaber_network.I, Kiefhaber_network.N, Kiefhaber_network.A], 
label=["I" "N" "A"], xlabel="Time (s)", ylabel="Concentration (M)")

# Compare MTK model and Reaction Network model
equations(Kiefhaber_network)
equations(convert(ODESystem, Kiefhaber_network))
equations(kiefhaber_model)
