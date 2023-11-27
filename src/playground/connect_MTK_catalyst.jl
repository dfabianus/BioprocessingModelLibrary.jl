# https://discourse.julialang.org/t/combining-catalyst-dsl-models-and-modelingtoolkit-explicit-ode-models/54418/13

using OrdinaryDiffEq
using ModelingToolkit
using Catalyst
using Plots

@parameters t k1
@variables A(t)
@derivatives D' ~ t
eqns = [
	D(A) ~ -k1*A
]

# Convert the system of equations to an ODESystem object
mdl1 = ODESystem(eqns, name=:mdl1)

mdl = @reaction_network begin
    @variables A(t) # This need to be explicitly declared as a variable
	k2*A,  0 --> B
end #k2 A

# Convert the system of reactions to and ODESystem object
tmp  = convert(ODESystem, mdl)
#mdl2 = ODESystem(tmp.eqs, pins = [A], name=:mdl2)
mdl2 = ODESystem(tmp.eqs, name=:mdl2)

# Connect the two systems
connections = [mdl2.A ~ mdl1.A]

# Create the combined model
@named combined = ODESystem(connections,t,[],[],systems=[mdl1,mdl2])
combined_simp = structural_simplify(combined)
# Set up parameters and initial conditions
u0 = [
	mdl1.A => 10,
	mdl2.B => 0
]

p = [
	mdl1.k1 => 2.0,
	mdl2.k2 => 1.0,
	#mdl2.A => 0
]

tspan = (0.0,30.0)

#prob = ODEProblem(alias_elimination(combined), u0, tspan, p)
prob = ODEProblem(alias_elimination(combined_simp), u0, tspan, p)
sol = solve(prob,Rodas5())

plot(sol)