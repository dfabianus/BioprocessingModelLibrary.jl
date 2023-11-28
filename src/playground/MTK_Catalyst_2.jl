using OrdinaryDiffEq
using ModelingToolkit
using Catalyst
using Plots

@parameters k1
@variables t A(t)
D = Differential(t)
eqns = [
	D(A) ~ -k1*A
    (@reaction k2*$A, 0 --> B)
]

@named mdl = ReactionSystem(eqns, t)
mdl = complete(mdl)


# Set up parameters and initial conditions
u0 = [
	mdl.A => 10,
	mdl.B => 0
]

p = [
	mdl.k1 => 2.0,
	mdl.k2 => 1.0,
]

tspan = (0.0,30.0)

#prob = ODEProblem(alias_elimination(combined), u0, tspan, p)
prob = ODEProblem(mdl, u0, tspan, p)
sol = solve(prob,Rodas5())

plot(sol)