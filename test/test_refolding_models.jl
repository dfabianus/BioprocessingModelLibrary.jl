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
    display(plot(sol_1, idxs=[sys.I, sys.N, sys.A], label=["I" "N" "A"], title="Testset Kiefhaber (Test 1: Standard simulation)", xlabel="Time (h)", ylabel="Concentration (M)"))

    kiefhaber_prob_2 = ODEProblem(sys, 
    [sys.I=>1.0±0.1, sys.N=>0.00±0.0, sys.A=>0.0±0.0], (0.0, 10.0),
    [sys.k_n=>0.2, sys.k_a=>0.1])
    sol_2 = solve(kiefhaber_prob_2, Tsit5(), saveat=0.1)
    @test sol_2.retcode == Success
    ts = range(tspan..., length=100)
    p = plot(ts, sol_2(ts, idxs=1).u, label="I(t)", title="Testset Kiefhaber (Test 2: MonteCarloMeasurements)", xlabel="Time (h)", ylabel="Concentration (M)")
    plot!(p, ts, sol_2(ts, idxs=2).u, label="N(t)")
    plot!(p, ts, sol_2(ts, idxs=3).u, label="A(t)")
    display(p)

    kiefhaber_prob_3 = ODEProblem(sys, 
    [sys.I=>1.0±0.1, sys.N=>0.00±0.0, sys.A=>0.0±0.0], (0.0, 10.0),
    [sys.k_n=>0.2±0.03, sys.k_a=>0.1±0.015])
    sol_3 = solve(kiefhaber_prob_3, Tsit5(), saveat=0.1)
    @test sol_3.retcode == Success
    ts = range(tspan..., length=100)
    p = plot(ts, sol_3(ts, idxs=1).u, label="I(t)", title="Testset Kiefhaber (Test 3: MonteCarloMeasurements)", xlabel="Time (h)", ylabel="Concentration (M)")
    plot!(p, ts, sol_3(ts, idxs=2).u, label="N(t)")
    plot!(p, ts, sol_3(ts, idxs=3).u, label="A(t)")
    display(p)
end

@testset "test_Kiefhaber_B" begin
    @named kiefhaber_model = Kiefhaber_B()
    sys = structural_simplify(kiefhaber_model)

    prob = ODEProblem(sys, 
    [sys.I=>1.0±0.1, sys.N=>0.00±0.0, sys.A=>0.0±0.0], (0.0, 10.0),
    [sys.a_n => 1.3343±0.0, sys.a_a => 12.0465±0.0, sys.b_n => -8.6824±0.7182, sys.b_a => -16.7869±2.5716, sys.D => 0.1])
    sol = solve(prob, Tsit5(), saveat=0.1)
    @test sol.retcode == Success
    ts = range(tspan..., length=100)
    p = plot(ts, sol(ts, idxs=1).u, label="I(t)", title="Testset Kiefhaber (Test 4: Rate dependency on D)", xlabel="Time (h)", ylabel="Concentration (M)")
    plot!(p, ts, sol(ts, idxs=2).u, label="N(t)")
    plot!(p, ts, sol(ts, idxs=3).u, label="A(t)")
    display(p)
end

# Reaction network models from Catalyst.jl are based on molar concentration, not mass concentrations.
@testset "test_Kiefhaber_network" begin
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
    xlabel="Time (h)", ylabel="Concentration (M)", title="Testset Kiefhaber Network (Test 1)")
    plot!(p, ts, sol(ts, idxs=3).u, label = "A(t)")
    display(p)
end


# Compare MTK model and Reaction Network model
# equations(kiefhaber_network)
# parameters(kiefhaber_network)
# states(kiefhaber_network)
# equations(convert(ODESystem, kiefhaber_network))
# equations(kiefhaber_model)
