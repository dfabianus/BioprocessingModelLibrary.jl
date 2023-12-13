# https://discourse.julialang.org/t/2d-arrays-with-modelingtoolkit/107448/6
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra

G_EARTH     = Float64[0.0, 0.0, -9.81]          # gravitational acceleration [m/s²]
L0::Float64 = -10.0                             # initial spring length      [m]
V0::Float64 = 4                                 # initial velocity           [m/s]

# model, Z component upwards
function springmass(;name)
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
    return ODESystem(eqs, t; name=name)
end

n=5
sys_array = []
for i in range(1, n)
    push!(sys_array, springmass(name=Symbol("springmass$i")))
end
for i in sys
    # Connect the springs in a chain
end



# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
# for l < l_0) and n tether segments. 
using ModelingToolkit, OrdinaryDiffEq, Plots, LinearAlgebra

G_EARTH     = Float64[0.0, 0.0, -9.81]          # gravitational acceleration     [m/s²]
L0::Float64 = 10.0                              # initial segment length            [m]
V0::Float64 = 2                                 # initial velocity of lowest mass [m/s]
segments::Int64 = 4                             # number of tether segments         [-]
α0 = π/8                                        # initial angle                   [rad]
POS0 = zeros(3, segments+1)
VEL0 = zeros(3, segments+1)
ACC0 = zeros(3, segments+1)
SEGMENTS0 = zeros(3, segments) 
UNIT_VECTORS0 = zeros(3, segments)
for i in 1:segments+1
    l0 = -(i-1)*L0
    v0 = (i-1)*V0/segments
    POS0[:, i] .= [sin(α0) * l0, 0, cos(α0) * l0]
    VEL0[:, i] .= [sin(α0) * v0, 0, cos(α0) * v0]
end
for i in 2:segments+1
    ACC0[:, i] .= G_EARTH
end
for i in 1:segments
    UNIT_VECTORS0[:, i] .= [0, 0, 1.0]
    SEGMENTS0[:, i] .= POS0[:, i+1] - POS0[:, i]
end

# model, Z component upwards
@parameters mass=1.0 c_spring0=50.0 damping=0.5 l_seg=L0
@variables t 
@variables pos(t)[1:3, 1:segments+1]  = POS0
@variables vel(t)[1:3, 1:segments+1]  = VEL0
@variables acc(t)[1:3, 1:segments+1]  = ACC0
@variables segment(t)[1:3, 1:segments]  = SEGMENTS0
@variables unit_vector(t)[1:3, 1:segments]  = UNIT_VECTORS0
@variables norm1(t)[1:segments] = l_seg * ones(segments)
@variables rel_vel(t)[1:3, 1:segments]  = zeros(3, segments)
@variables spring_vel(t)[1:segments] = zeros(segments)
@variables c_spring(t)[1:segments] = c_spring0 * ones(segments)
@variables spring_force(t)[1:3, 1:segments] = zeros(3, segments)
@variables total_force(t)[1:3, 1:segments] = zeros(3, segments)
D = Differential(t)

eqs1 = vcat(D.(pos) ~ vel,
            D.(vel) ~ acc)
eqs2 = []
for i in segments:-1:1
    global eqs2
    eqs2 = vcat(eqs2, segment[:, i] ~ pos[:, i+1] - pos[:, i])
    eqs2 = vcat(eqs2, norm1[i] ~ norm(segment[:, i]))
    eqs2 = vcat(eqs2, unit_vector[:, i] ~ -segment[:, i]/norm1[i])
    eqs2 = vcat(eqs2, rel_vel[:, i] ~ vel[:, i+1] - vel[:, i])
    eqs2 = vcat(eqs2, spring_vel[i] ~ -unit_vector[:, i] ⋅ rel_vel[:, i])
    eqs2 = vcat(eqs2, c_spring[i] ~ c_spring0 * (norm1[i] > l_seg))
    eqs2 = vcat(eqs2, spring_force[:, i] ~ (c_spring[i] * (norm1[i] - l_seg) + damping * spring_vel[i]) * unit_vector[:, i])
    if i == segments
        eqs2 = vcat(eqs2, total_force[:, i] ~ spring_force[:, i])
    else
        eqs2 = vcat(eqs2, total_force[:, i] ~ spring_force[:, i]- spring_force[:, i+1])
    end
    eqs2 = vcat(eqs2, acc[:, i+1] .~ G_EARTH + total_force[:, i] / mass)
end
eqs2 = vcat(eqs2, acc[:, 1] .~ zeros(3))
eqs = vcat(eqs1..., eqs2)
     
@named sys = ODESystem(eqs, t)
simple_sys = structural_simplify(sys)

duration = 10.0
dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration

u0 = Float64[]
for i in 1:segments+1
    global u0
    u0 = vcat(u0, POS0[:, i], VEL0[:, i])
end

prob = ODEProblem(simple_sys, u0, tspan)
@time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)

function plot2d(sol, reltime, segments)
    index = Int64(round(reltime*50+1))
    x = Float64[]
    z = Float64[]
    for particle in 1:segments+1
        push!(x, (sol(sol.t, idxs=pos[1, particle]))[index])
        push!(z, (sol(sol.t, idxs=pos[3, particle]))[index])
    end
    x_max = maximum(x)
    z_max = maximum(z)
    plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false)
    annotate!(15, z_max-3.0, "t=$(round(reltime,digits=1)) s")
    plot!(x, z, seriestype = :scatter) 
    ylims!((-segments*10-10, 0.5))
    xlims!((-segments*5, segments*5))
end

dt = 0.04
for time in 0:dt:10
    display(plot2d(sol, time, segments))
    sleep(0.25*dt)
end
nothing



# Example two: Falling mass, attached to non-linear spring without compression stiffness,
# initially moving upwards with 4 m/s, using a callback to precisely calculate the discontinuities
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra

G_EARTH  = Float64[0.0, 0.0, -9.81]             # gravitational acceleration [m/s²]
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

@named sys = ODESystem(eqs, t; continuous_events = [norm(pos) ~ abs(L0)])
simple_sys = structural_simplify(sys)

duration = 10.0
dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration
# initial state
u0 = Dict(pos=>[0,0,L0], vel=>[0,0,V0])

# function condition(u, t, integrator) # Event when condition(u,t,integrator) == 0
#     norm(u[1:3]) - abs(L0)
# end
# function affect!(integrator)
#     println(integrator.t)
# end
cb = ContinuousCallback(condition, affect!)

prob = ODEProblem(simple_sys, u0, tspan)
sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts, callback = cb)

X = sol.t
POS_Z = sol(X, idxs=pos[3])
VEL_Z = sol(X, idxs=vel[3])

lns1 = plot(X, POS_Z, color="green", label="pos_z")
xlabel("time [s]")
ylabel("pos_z [m]")
lns2 = plot(X, L0.+0.005 .* sol[c_spring], color="grey", label="c_spring")
grid(true)
twinx()
ylabel("vel_z [m/s]") 
lns3 = plot(X, VEL_Z, color="red", label="vel_z")
lns = vcat(lns1, lns2, lns3)
labs = [l.get_label() for l in lns]
legend(lns, labs) 

nothing