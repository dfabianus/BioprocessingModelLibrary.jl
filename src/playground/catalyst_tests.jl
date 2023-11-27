@parameters t
Dt = Differential(t)

function TEST_Enzymkinetik()
    TEST_Enzymkinetik = @reaction_network begin
        (k1,k1r), E + S <--> ES
        k2, ES --> E + P
    end
    return TEST_Enzymkinetik
end

model = TEST_Enzymkinetik()
equations(model)
equations(convert(ODESystem, model))
p = (:k1 => 1000.0, :k1r => 0.00001, :k2 => 1.0)
tspan = (0.,20.)
u0 = [:E => 1.0, :S => 10.0, :ES => 0.0, :P => 0.0]
oprob = ODEProblem(model, u0, tspan, p)
osol  = solve(oprob, Tsit5())
plot(osol)
km = 0.00001/1000.0
plot(osol, idxs = [:E, :S, :ES, :P], title = "Enzymkinetik", xlabel = "Zeit", ylabel = "Konzentration", label = ["E" "S" "ES" "P"])
plot(osol, idxs = (:S))


# It works
function FLUMO3()
    FLUMO3 = @reaction_network begin
        k_a, 2I --> A
        k_ac, I + C --> A
        k_ic, I + C --> IC
        k_nc, IC --> NC
        k_n, I --> N
        k_cn, N + C --> NC
    end
    return FLUMO3
end

@mtkmodel FLUMO_D begin
    @parameters begin
        c_Din
    end
    @variables begin
        F_in(t)
        D(t)
    end
    @equations begin
        Dt(D) ~ F_in * c_Din
        F_in ~ 0.1
    end
end

@named model1 = FLUMO_D()
model2 = convert(ODESystem, FLUMO3())
@named combined = extend(model1, model2)
states(combined)
equations(combined)
parameters(combined)
sys = structural_simplify(combined)
p = (sys.k_a => 1.0, sys.k_ac => 1.0, sys.k_ic => 1.0, sys.k_nc => 1.0, sys.k_n => 1.0, sys.k_cn => 1.0, sys.c_Din => 1.0)
tspan = (0.,20.)
u0 = [sys.D => 1.0, sys.I => 1.0, sys.A => 1.0, sys.C => 1.0, sys.IC => 1.0, sys.N => 1.0, sys.NC => 1.0]
oprob = ODEProblem(sys, u0, tspan, p)
osol  = solve(oprob, Tsit5())
plot(osol)

# test stoichiomatric matrix and complex stoichiomatric matrix
function TEST_SMATRIX()
    rn = @reaction_network begin
        k1, 6C + 12H + 6O --> Glc
        k2, 6C + 12H + 6O --> X
    end
    return rn
end
TEST_SMATRIX()
smatrix = netstoichmat(TEST_SMATRIX())
smatrix = complexstoichmat(TEST_SMATRIX())
species(TEST_SMATRIX())
reactioncomplexes(TEST_SMATRIX())