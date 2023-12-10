function Pablymueller(induction_time)
    @parameters a_n b_n a_a b_a a_cn a_ic a_nc p1 p2 p3 p4 p5
    @variables begin
        P(t)
        V(t)
        D(t)
        I(t)
        N(t)
        A(t)
        NC(t)
        IC(t)
        C(t)
        k_a(t)
        k_n(t)
        k_cn(t)
        k_ic(t)
        k_nc(t)
        cI(t)
        cN(t)
        cA(t)
        cD(t)
        cP(t)
        cNC(t)
        cIC(t)
        cC(t)
        AEW(t)#
        F(t)
    end

    eqns = [
        Dt(V) ~ 0
        Dt(D) ~ 0
        P ~ I + N + NC + IC # +2*A
        k_a ~ maximum([0, a_a * (1 + D) ^ b_a])
        k_n ~ maximum([0, a_n * (1 + D) ^ b_n])
        k_cn ~ maximum([0, a_cn])
        k_ic ~ maximum([0, a_ic])
        k_nc ~ maximum([0, a_nc])
        cI  ~ I / V
        cN  ~ N / V
        cA  ~ A / V
        cD  ~ D / V
        cP  ~ P / V
        cNC ~ NC / V
        cIC ~ IC / V
        cC  ~ C / V
        F ~ p1 * cP / (p2 + cP)
        AEW ~ -(abs(p3) * ((cN+cA)/cP) / (abs(p4) + ((cN+cA)/cP)) + abs(p5)*((cN+cA)/cP))
        
        (@reaction $k_a,  2*I --> A)
        (@reaction $k_n,  I --> N)
        (@reaction $k_cn, N + C --> NC)
        (@reaction $k_ic, I + C --> IC)
        (@reaction $k_nc, IC --> NC)
    ]

    mds = D .~ D .+ mD_pulses
    mps = I .~ I .+ mP_pulses
    mcs = C .~ C .+ mC_pulses
    Vs =  V .~ V .+ V_pulses

    @named mdl = ReactionSystem(eqns, t; discrete_events = vcat(
        [[a] => [b] for (a,b) in zip(pulse_times, mds)],
        [[a] => [b] for (a,b) in zip(pulse_times, mps)],
        [[a] => [b] for (a,b) in zip(pulse_times, mcs)],
        [[a] => [b] for (a,b) in zip(pulse_times, Vs)],
    ))

    return structural_simplify(convert(ODESystem, mdl))
end