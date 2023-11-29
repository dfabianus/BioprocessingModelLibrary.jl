using OrdinaryDiffEq
using ModelingToolkit
using Catalyst
using Plots

@parameters k1
@variables t A(t) A2(t)
D = Differential(t)
eqns = [
	D(A) ~ -k1*A
	A2 ~ A * 2
    (@reaction k2*$A2, 0 --> B)
]

@named mdl = ReactionSystem(eqns, t)
mdl_simp = structural_simplify(convert(ODESystem, mdl))
#complete(mdl)

# Set up parameters and initial conditions
u0 = [
	mdl_simp.A => 10,
	mdl_simp.B => 0
]

p = [
	mdl_simp.k1 => 2.0,
	mdl_simp.k2 => 1.0,
]

tspan = (0.0,30.0)

prob = ODEProblem(mdl_simp, u0, tspan, p)
sol = solve(prob,Rodas5())

plot(sol)