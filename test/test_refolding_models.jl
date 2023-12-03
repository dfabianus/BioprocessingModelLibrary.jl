using Revise
using ModelingToolkit
using OrdinaryDiffEq
using OrdinaryDiffEq: ReturnCode.Success
using BioprocessingModelLibrary.Refolding
using MonteCarloMeasurements
using Plots
using Test

@testset "test_Kiefhaber" begin
    tspan = (0.0, 10.0)
    @named kiefhaber_model = Kiefhaber()
    sys = structural_simplify(kiefhaber_model)

    kiefhaber_prob_1 = ODEProblem(sys, 
    [sys.I=>1.0, sys.N=>0.0, sys.A=>0.0], (0.0, 10.0),
    [sys.k_n=>0.2, sys.k_a=>0.1])
    sol_1 = solve(kiefhaber_prob_1, Tsit5(), saveat=0.1)
    @test sol_1.retcode == Success
    display(plot(sol_1, idxs=[sys.I, sys.N, sys.A], label=["I" "N" "A"], title="Testset Kiefhaber (Test 1: Standard simulation)", xlabel="Time (h)", ylabel="Concentration (g/L)"))

    kiefhaber_prob_2 = ODEProblem(sys, 
    [sys.I=>1.0±0.1, sys.N=>0.00±0.0, sys.A=>0.0±0.0], (0.0, 10.0),
    [sys.k_n=>0.2, sys.k_a=>0.1])
    sol_2 = solve(kiefhaber_prob_2, Tsit5(), saveat=0.1)
    @test sol_2.retcode == Success
    ts = range(tspan..., length=100)
    p = plot(ts, sol_2(ts, idxs=1).u, label="I(t)", title="Testset Kiefhaber (Test 2: MonteCarloMeasurements)", xlabel="Time (h)", ylabel="Concentration (g/L)")
    plot!(p, ts, sol_2(ts, idxs=2).u, label="N(t)")
    plot!(p, ts, sol_2(ts, idxs=3).u, label="A(t)")
    display(p)

    kiefhaber_prob_3 = ODEProblem(sys, 
    [sys.I=>1.0±0.1, sys.N=>0.00±0.0, sys.A=>0.0±0.0], (0.0, 10.0),
    [sys.k_n=>0.2±0.03, sys.k_a=>0.1±0.015])
    sol_3 = solve(kiefhaber_prob_3, Tsit5(), saveat=0.1)
    @test sol_3.retcode == Success
    ts = range(tspan..., length=100)
    p = plot(ts, sol_3(ts, idxs=1).u, label="I(t)", title="Testset Kiefhaber (Test 3: MonteCarloMeasurements)", xlabel="Time (h)", ylabel="Concentration (g/L)")
    plot!(p, ts, sol_3(ts, idxs=2).u, label="N(t)")
    plot!(p, ts, sol_3(ts, idxs=3).u, label="A(t)")
    display(p)
end

@testset "test_Kiefhaber_B" begin
    tspan = (0.0, 10.0)
    @named kiefhaber_model = Kiefhaber_B()
    sys = structural_simplify(kiefhaber_model)
    prob = ODEProblem(sys, 
    [sys.I=>1.0±0.1, sys.N=>0.00±0.0, sys.A=>0.0±0.0], (0.0, 10.0),
    [sys.a_n => 1.3343±0.0, sys.a_a => 12.0465±0.0, sys.b_n => -8.6824±0.7182, sys.b_a => -16.7869±2.5716, sys.D => 0.1])
    sol = solve(prob, Tsit5(), saveat=0.1)
    @test sol.retcode == Success
    ts = range(tspan..., length=100)
    p = plot(ts, sol(ts, idxs=1).u, label="I(t)", title="Testset Kiefhaber (Test 4: Rate dependency on D)", xlabel="Time (h)", ylabel="Concentration (g/L)")
    plot!(p, ts, sol(ts, idxs=2).u, label="N(t)")
    plot!(p, ts, sol(ts, idxs=3).u, label="A(t)")
    display(p)
end

# Reaction network models from Catalyst.jl are based on molar concentration, not mass concentrations.
@testset "test_Kiefhaber_network" begin
    tspan = (0.0, 10.0)
    kiefhaber_network = Kiefhaber_network()
    prob = ODEProblem(kiefhaber_network, 
        [:I => 1.0±0.1, :N => 0.0±0.0, :A => 0.0±0.0], 
        tspan, 
        (:a_n => 1.3343±0.0, :a_a => 12.0465±0.0, :b_n => -8.6824±0.7182, :b_a => -16.7869±2.5716, :D => 0.1))
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    @test sol.retcode == Success
    ts = range(tspan..., length=100)
    p = plot(ts, sol(ts, idxs=1).u, label = "I(t)")
    plot!(p, ts, sol(ts, idxs=2).u, label = "N(t)",
    xlabel="Time (h)", ylabel="Concentration (mol/L)", title="Testset Kiefhaber Network (Test 1)")
    plot!(p, ts, sol(ts, idxs=3).u, label = "A(t)")
    display(p)
end


@testset "test_INAM" begin
    tspan = (0.0, 10.0)
    @named sys = INAM()
    sys_simp = structural_simplify(sys)
    prob = ODEProblem(sys_simp, 
    [sys_simp.I=>1.0±0.1, sys_simp.N=>0.00±0.0, sys_simp.A=>0.0±0.0, sys_simp.M=>0.0±0.0], (0.0, 10.0),
    [sys_simp.a_n => 1.3343±0.0, sys_simp.a_a => 12.0465±0.0, sys_simp.a_m => -1.0±0,
    sys_simp.b_n => -8.6824±0.7182, sys_simp.b_a => -16.7869±2.5716,  sys_simp.b_m => 1.0±0,
    sys_simp.D => 0.1])
    sol = solve(prob, Tsit5(), saveat=0.1)
    @test sol.retcode == Success
    ts = range(tspan..., length=100)
    p = plot(ts, sol(ts, idxs=1).u, label="I(t)", title="Testset INAM (Test 1: The naive implementation)", xlabel="Time (h)", ylabel="Concentration (g/L)")
    plot!(p, ts, sol(ts, idxs=2).u, label="N(t)")
    plot!(p, ts, sol(ts, idxs=3).u, label="A(t)")
    plot!(p, ts, sol(ts, idxs=4).u, label="M(t)")
    display(p)
end

@testset "test_INAM_network" begin
    tspan = (0.0, 10.0)
    rn = INAM_network()
    prob = ODEProblem(rn, 
        [:I => 1.0±0.1, :N => 0.0±0.0, :A => 0.0±0.0, :M => 0.0±0.0], 
        tspan, 
        (:a_n => 1.3343±0.0, :a_a => 12.0465±0.0, :a_m => 12.0465±0.0,
        :b_n => -8.6824±0.7182, :b_a => -16.7869±2.5716, :b_m => -16.7869±2.5716,
        :D => 0.1)
    )
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    @test sol.retcode == Success
    ts = range(tspan..., length=100)
    p = plot(ts, sol(ts, idxs=1).u, label = "I(t)")
    plot!(p, ts, sol(ts, idxs=2).u, label = "N(t)",
    xlabel="Time (h)", ylabel="Concentration (mol/L)", title="Testset Kiefhaber Network (Test 1)")
    plot!(p, ts, sol(ts, idxs=3).u, label = "A(t)")
    plot!(p, ts, sol(ts, idxs=4).u, label = "M(t)")
    display(p)
end

@testset "test_MTK_Catalyst_Connections" begin
    sys = FLUMO_COMBINED_3()
    p = (sys.k_a => 1.0, sys.k_ac => 1.0, sys.k_ic => 1.0, sys.k_nc => 1.0, 
    sys.k_n => 1.0, sys.k_cn => 1.0, sys.c_Din => 1.0)
    tspan = (0.,20.)
    u0 = [sys.D => 1.0, sys.I => 1.0, sys.A => 1.0, sys.C => 1.0, sys.IC => 1.0, 
    sys.N => 1.0, sys.NC => 1.0]
    oprob = ODEProblem(sys, u0, tspan, p)
    osol  = solve(oprob, Tsit5())
    @test osol.retcode == Success
    
    sys = FLUMO_COMBINED_4()
    p = (sys.reactor.c_Din => 5, sys.reactions_ODE.k_a => 1.0, sys.reactions_ODE.k_c => 1.0, 
    sys.reactions_ODE.a_n => 1.0, sys.reactions_ODE.b_n => 1.0)
    tspan = (0.,20.)
    u0 = [sys.reactor.D => 1.0, sys.reactions_ODE.I => 1.0, sys.reactions_ODE.A => 1.0, 
    sys.reactions_ODE.C => 1.0, sys.reactions_ODE.IC => 1.0, sys.reactions_ODE.N => 1.0]
    oprob = ODEProblem(sys, u0, tspan, p)
    osol  = solve(oprob, Tsit5())
    @test osol.retcode == Success

    sys = FLUMO_COMBINED_5()
    p = (sys.reactor.c_Din => 5, sys.kinetics.a_n => 1.0, sys.kinetics.b_n => 1.0, 
        sys.kinetics.a_a => 1.0, sys.kinetics.a_ac => 1.0, sys.kinetics.a_ic => 1.0,
        sys.kinetics.a_nc => 1.0, sys.kinetics.a_cn => 1.0
    )
    u0 = [sys.reactor.D => 1.0, sys.reactions_ODE.I => 1.0, sys.reactions_ODE.A => 0.0, 
        sys.reactions_ODE.C => 0.5, sys.reactions_ODE.IC => 0.0, sys.reactions_ODE.NC => 0.0,
        sys.reactions_ODE.N => 0.0
    ]
    tspan = (0.,5.)
    oprob = ODEProblem(sys, u0, tspan, p)
    osol  = solve(oprob, Tsit5())
    @test osol.retcode == Success
    p = plot(osol, title = "Testset MTK Catalyst Connections", xlabel="Time (h)", ylabel="Concentration (mol/L)")
    display(p)
end

@testset "test_FLUMO_pulse_batch_processing" begin
    pulse_times = [1.0, 2.0, 3.0, 4.0, 5.0]
    mP_pulses = [0.1, 0.2, 0.2, 0.5, 0.1]
    mD_pulses = [0.1, 0.1, 0.1, 0.1, 0.1]
    mC_pulses = [0.4, 0.2, 0.2, 0.1, 0.3]
    V_pulses = [0.2, 0.2, 0.2, 0.2, 0.2]
    sys = FLUMO_COMBINED_6(pulse_times, mP_pulses, mD_pulses, mC_pulses, V_pulses)
    p = (sys.reactor.c_Din => 5, sys.kinetics.a_n => 1.0, sys.kinetics.b_n => 1.0, 
        sys.kinetics.a_a => 1.0, sys.kinetics.a_ac => 1.0, sys.kinetics.a_ic => 1.0,
        sys.kinetics.a_nc => 1.0, sys.kinetics.a_cn => 1.0
    )
    u0 = [sys.reactor.D => 1.0, sys.reactions_ODE.I => 1.0, sys.reactions_ODE.A => 0.0, 
        sys.reactions_ODE.C => 0.5, sys.reactions_ODE.IC => 0.0, sys.reactions_ODE.NC => 0.0,
        sys.reactions_ODE.N => 0.0, sys.reactor.V => 1.0
    ]
    tspan = (0.,5.)
    oprob = ODEProblem(sys, u0, tspan, p)
    osol  = solve(oprob, Tsit5())
    @test osol.retcode == Success
    p = plot(osol, title = "Testset MTK Catalyst Connections", xlabel="Time (h)", ylabel="Concentration (mol/L)")
    display(p)
end

@testset "test_FLUMO_pulse_batch_processing" begin
    pulse_times = [1.0, 2.0, 3.0, 4.0, 5.0]
    mP_pulses = [0.1, 0.2, 0.2, 0.5, 0.1]
    mD_pulses = [0.1, 0.1, 0.1, 0.1, 0.1]
    mC_pulses = [0.4, 0.2, 0.2, 0.1, 0.3]
    V_pulses = [0.2, 0.2, 0.2, 0.2, 0.2]
    sys = FLUMO_COMBINED_6(pulse_times, mP_pulses, mD_pulses, mC_pulses, V_pulses)
    p = (sys.reactor.c_Din => 5, sys.kinetics.a_n => 1.0, sys.kinetics.b_n => 1.0, 
        sys.kinetics.a_a => 1.0, sys.kinetics.b_a => 1.0, sys.kinetics.a_ac => 1.0, sys.kinetics.a_ic => 1.0,
        sys.kinetics.a_nc => 1.0, sys.kinetics.a_cn => 1.0
    )
    u0 = [sys.reactor.D => 1.0±0.1, sys.reactions_ODE.I => 1.0±0.1, sys.reactions_ODE.A => 0.0, 
        sys.reactions_ODE.C => 0.5±0.05, sys.reactions_ODE.IC => 0.0, sys.reactions_ODE.NC => 0.0,
        sys.reactions_ODE.N => 0.0, sys.reactor.V => 1.0±0.05
    ]
    tspan = (0.,5.)
    oprob = ODEProblem(sys, u0, tspan, p)
    osol  = solve(oprob, Tsit5())
    @test osol.retcode == Success
    ts = range(tspan..., length=1000)
    p = plot(ts, osol(ts, idxs=sys.reactions_ODE.I).u, label = "I(t)")
    p = plot!(ts, osol(ts, idxs=sys.reactions_ODE.A).u, label = "A(t)")
    display(p)
end

@testset "test_FLUMO_MECH_pulse_cofactor" begin
    pulse_times = [1.0, 2.0, 3.0, 4.0, 5.0]
    mP_pulses = [0.1, 0.2, 0.2, 0.5, 0.1]
    mD_pulses = [0.1, 0.1, 0.1, 0.1, 0.1]
    V_pulses = [0.2, 0.2, 0.2, 0.2, 0.2]
    sys = FLUMO_MECH(pulse_times, mP_pulses, mD_pulses, zeros(5), V_pulses)
    p = (sys.a_n => 1.0, sys.b_n => 1.0, 
        sys.a_a => 1.0, sys.b_a => 1.0, sys.a_cn => 1.0, sys.a_ic => 1.0, sys.a_nc => 1.0,
        sys.p1 => 1.0, sys.p2 => 1.0, sys.p3 => 1.0, sys.p4 => 1.0, sys.p5 => 1.0
    )
    u0 = [sys.D => 1.0±0.1, sys.I => 1.0±0.1, sys.A => 0.0, 
        sys.N => 0.0, sys.V => 1.0±0.05, sys.C => 0.1±0.01, sys.IC => 0, sys.NC => 0
    ]
    tspan = (0.,6.)
    oprob = ODEProblem(sys, u0, tspan, p)
    osol  = solve(oprob, Tsit5())
    @test osol.retcode == Success
    ts = range(tspan..., length=1000)
    p = plot(ts, osol(ts, idxs=sys.I).u, label = "I(t)")
    p = plot!(ts, osol(ts, idxs=sys.A).u, label = "A(t)")
    p = plot!(ts, osol(ts, idxs=sys.IC).u, label = "IC(t)")
    p = plot!(ts, osol(ts, idxs=sys.NC).u, label = "NC(t)")
    p = plot!(ts, osol(ts, idxs=sys.cP+2*sys.cA).u, label = "cP(t)")
    display(p)
    p2 = plot(ts, osol(ts, idxs=sys.AEW).u, label = "AEW(t)")
    display(p2)
end


    @named sys = FLUMO_SOFTSENSOR()
    sys_simp = structural_simplify(sys)
    p = (sys_simp.f_IAEW => 1.0, sys_simp.p1 => 1.0, 
        sys_simp.p2 => 1.0, 
    )
    u0 = [sys_simp.I => 1.0, sys_simp.N => 0.0,
    ]
    tspan = (0.,6.)
    oprob = ODEProblem(sys_simp, u0, tspan, p)
    osol  = solve(oprob)
    @test osol.retcode == Success
    ts = range(tspan..., length=1000)
    p = plot(ts, osol(ts, idxs=sys_simp.I).u, label = "I(t)")
    p = plot(ts, osol(ts, idxs=sys_simp.A).u, label = "A(t)", ylim=(0,150))
    p = plot(ts, osol(ts, idxs=sys_simp.N).u, label = "N(t)")