@mtkmodel FLUMO_ReactorDynamics1 begin
    @parameters begin
        c_Din
    end
    @variables begin
        F_in(t)
        D(t)
    end
    @equations begin
        Dt(D) ~ F_in * c_Din
        F_in ~ 0.0
    end
end

@mtkmodel FLUMO_ReactorDynamics2 begin
    @parameters begin
        c_Din
    end
    @variables begin
        F_in(t)
        D(t)
        V(t)
    end
    @equations begin
        Dt(V) ~ F_in
        Dt(D) ~ F_in * c_Din
        F_in ~ 0.0
    end
end

@mtkmodel FLUMO_ReactorDynamics3 begin
    @parameters begin
        c_Din
    end
    @variables begin
        F_in(t)
        V(t)
        I(t)
        N(t)
        A(t)
        D(t)
        P(t)
        cI(t)
        cN(t)
        cA(t)
        cD(t)
        cP(t)
    end
    @equations begin
        Dt(V) ~ F_in
        Dt(D) ~ F_in * c_Din
        P ~ I + N
        F_in ~ 0.0
        cI  ~ I / V
        cN  ~ N / V
        cA  ~ A / V
        cD  ~ D / V
        cP  ~ P / V
    end
end

@mtkmodel FLUMO_Kinetics1 begin
    @parameters begin
        a_n
        b_n
        a_a
        a_ac 
        a_ic
        a_nc
        a_cn 
    end
    @variables begin
        k_n(t)
        k_a(t)
        k_ac(t)
        k_ic(t)
        k_nc(t)
        k_cn(t)
        D(t)
    end
    @equations begin
        k_n ~ a_n * (1 + D) ^ b_n
        k_a ~ a_a
        k_ac ~ a_ac
        k_ic ~ a_ic
        k_nc ~ a_nc
        k_cn ~ a_cn
    end
end

function FLUMO_Reactions1()
    FLUMO1 = @reaction_network begin
        k_a, n*I --> A
        k_c, I + C --> IC
        k_n, IC --> N
    end
    return FLUMO1
end

function FLUMO_Reactions2()
    FLUMO2 = @reaction_network begin
        k_a, n*I --> A
        k_ac, I + C --> A
        k_nc, I + C --> IC
        k_n, IC --> N
    end
    return FLUMO2
end

function FLUMO_Reactions3()
    FLUMO3 = @reaction_network begin
        k_a, 2*I --> A
        k_ac, I + C --> A
        k_ic, I + C --> IC
        k_nc, IC --> NC
        k_n, I --> N
        k_cn, N + C --> NC
    end
    return FLUMO3
end

function FLUMO_COMBINED_3()
    @named reactor = FLUMO_ReactorDynamics1()
    reactions = convert(ODESystem, FLUMO_Reactions3())
    @named FLUMO_COMBINED_3 = extend(reactor, reactions)
    return structural_simplify(FLUMO_COMBINED_3)
end

function FLUMO_Reactions4()
    FLUMO4 = @reaction_network begin
        @variables D(t)
        k_a, 2*I --> A
        k_c, I + C --> IC
        a_n*(1+D)^b_n, IC --> N
    end
    return FLUMO4
end

FLUMO_Reactions5 = @reaction_network begin
    @variables k_a(t) k_ac(t) k_ic(t) k_nc(t) k_n(t) k_cn(t)
    k_a, 2*I --> A
    k_ac, I + C --> A
    k_ic, I + C --> IC
    k_nc, IC --> NC
    k_n, I --> N
    k_cn, N + C --> NC
end

function FLUMO_COMBINED_4()
    @named reactor = FLUMO_ReactorDynamics1()
    reactions = convert(ODESystem, FLUMO_Reactions4())
    @named reactions_ODE = ODESystem(reactions.eqs)
    connections = [reactor.D ~ reactions_ODE.D]
    @named combined = ODESystem(connections,t,[],[],systems=[reactor, reactions_ODE])
    return structural_simplify(combined)
end

function FLUMO_COMBINED_5()
    @named reactor = FLUMO_ReactorDynamics1()
    @named kinetics = FLUMO_Kinetics1()
    reactions = convert(ODESystem, FLUMO_Reactions5)
    @named reactions_ODE = ODESystem(reactions.eqs)
    connections = [
        reactor.D ~ kinetics.D,
        kinetics.k_a ~ reactions_ODE.k_a,
        kinetics.k_ac ~ reactions_ODE.k_ac,
        kinetics.k_ic ~ reactions_ODE.k_ic,
        kinetics.k_nc ~ reactions_ODE.k_nc,
        kinetics.k_n ~ reactions_ODE.k_n,
        kinetics.k_cn ~ reactions_ODE.k_cn
        ]
    @named combined = ODESystem(connections,t,[],[],systems=[reactor, kinetics, reactions_ODE])
    return structural_simplify(combined)
end

function FLUMO_COMBINED_6(pulse_times, mP_pulses, mD_pulses, mC_pulses, V_pulses)
    @named reactor = FLUMO_ReactorDynamics2()
    @named kinetics = FLUMO_Kinetics1()
    reactions = convert(ODESystem, FLUMO_Reactions5)
    @named reactions_ODE = ODESystem(reactions.eqs)
    times = pulse_times
    mds = reactor.D .~ reactor.D .+ mD_pulses
    mps = reactions_ODE.I .~ reactions_ODE.I .+ mP_pulses
    mcs = reactions_ODE.C .~ reactions_ODE.C .+ mC_pulses
    Vs =  reactor.V .~ reactor.V .+ V_pulses
    connections = [
        reactor.D ~ kinetics.D,
        kinetics.k_a ~ reactions_ODE.k_a,
        kinetics.k_ac ~ reactions_ODE.k_ac,
        kinetics.k_ic ~ reactions_ODE.k_ic,
        kinetics.k_nc ~ reactions_ODE.k_nc,
        kinetics.k_n ~ reactions_ODE.k_n,
        kinetics.k_cn ~ reactions_ODE.k_cn
    ]
    @named combined = ODESystem(connections,t,[],[],systems=[reactor, kinetics, reactions_ODE],
    discrete_events = vcat(
        [[a] => [b] for (a,b) in zip(times, mds)],
        [[a] => [b] for (a,b) in zip(times, mps)],
        [[a] => [b] for (a,b) in zip(times, mcs)],
        [[a] => [b] for (a,b) in zip(times, Vs)],
    ))
    return structural_simplify(combined)
end