@mtkmodel FLUMO_ReactorDynamics begin
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
    @named reactor = FLUMO_ReactorDynamics()
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

function FLUMO_COMBINED_4()
    @named reactor = FLUMO_ReactorDynamics()
    reactions = convert(ODESystem, FLUMO_Reactions4())
    @named reactions_ODE = ODESystem(reactions.eqs)
    connected = compose(ODESystem([reactor.D ~ reactions_ODE.D], t; name = :connect), reactor, reactions_ODE)
    return structural_simplify(connected)
end


