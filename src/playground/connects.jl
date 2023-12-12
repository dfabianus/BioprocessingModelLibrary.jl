using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra

G_EARTH     = Float64[0.0, 0.0, -9.81]          # gravitational acceleration [m/s²]
L0::Float64 = -10.0                             # initial spring length      [m]
V0::Float64 = 4                                 # initial velocity           [m/s]

# model, Z component upwards
@parameters mass=1.0 c_spring0=50.0 damping=0.5 l0=L0
@variables t pos(t)[1:3] = [0.0, 0.0,  L0]
@variables   vel(t)[1:3] = [0.0, 0.0,  V0] 
@variables   acc(t)[1:3] = G_EARTH
@variables unit_vector(t)[1:3]  = [0.0, 0.0, -sign(L0)]
@variables c_spring(t) = c_spring0
@variables spring_force(t)[1:3] = [0.0, 0.0, 0.0]
@variables force(t) = 0.0 norm1(t) = abs(l0) spring_vel(t) = 0.0
D = Differential(t)

eqs = vcat(D.(pos)      ~ vel,
           D.(vel)      ~ acc,
           norm1        ~ norm(pos),
           unit_vector  ~ -pos/norm1,         # direction from point mass to origin
           spring_vel   ~ -unit_vector ⋅ vel,
           c_spring     ~ c_spring0 * (norm1 > abs(l0)),
           spring_force ~ (c_spring * (norm1 - abs(l0)) + damping * spring_vel) * unit_vector,
           acc          ~ G_EARTH + spring_force/mass)