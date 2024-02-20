function ColiModel(induction_time, qSset, qSset2; name, beta2=0.)
    @parameters begin
        qSmax = 1.2
        Ks = 0.04
        tInd = induction_time
        Yxs = 0.5
        kd = 0.0
        cSR = 400.
        Ms = 30.
        Mx = 26.5
        DorS = 4.
        DorX = 4.114
        DorO2 = -4.
        qPmax   = 5.83E-2 #± 1.79E-4;
        qPIBmax = 2.14E-2 #± 2.78E-4;
        kSqS    = 14.5E-2 #± 1.13E-4;
        k       = 6.50 #± 1.06E-2;
        kIqP    = 6.03 #± 5.61E-2;
        kSqP    = 70.1E-2 #± 3.64E-3;
        N       = 1.496 #± 1.51E-3;
        kPqPIB  = 3.742E-2 #± 4.23E-4;
    end
    @variables begin
        V(t)
        Fin(t)
        Xv(t)
        Xd(t)
        S(t)
        CER(t)
        OUR(t)
        qS(t)
        qP(t)
        qPIB(t)
        qScap(t)
        β1(t) = 0.
        qSw(t) = 0.
        Ptot(t)
        PIB(t)
        Smet(t)
        Smetspec(t)
        STY(t)
    end

    eqns = [
        Dt(β1) ~ 0
        Dt(qSw) ~ 0
        Fin ~ qSw * Xv*V / cSR
        qScap ~ qSmax * exp(β1 * (t - tInd))
        qS ~ qScap * S / (Ks + S)
        Dt(V) ~ Fin 
        Dt(Xv) ~ Yxs * qS * Xv - kd * Xv - Fin/V * Xv
        Dt(Xd) ~ kd * Xd - Fin/V * Xd
        Dt(S) ~ Fin/V * cSR - Fin/V * S - qS * Xv
        Dt(Ptot) ~ IfElse.ifelse(t>induction_time, qP * Xv - kd * Ptot, 0.) 
        Dt(PIB) ~ IfElse.ifelse(t>induction_time, qPIB * Xv - kd * PIB, 0.)
        Dt(Smet) ~ IfElse.ifelse(t>induction_time, qS * Xv, 0.)
        STY ~ IfElse.ifelse(t>induction_time, qP*Xv, 0.)
        Smetspec ~ Smet / Xv
        CER ~ (1/Ms - Yxs/Mx)*qS*Xv
        OUR ~ (DorS/Ms - Yxs * DorX/Mx)*qS*Xv/(DorO2)
        qP ~ qPmax * qS / (kSqS + qS) * Smetspec / ((Smetspec^k /kIqP)+Smetspec+kSqP)
        qPIB ~ qPIBmax * qP^N / (qP^N + kPqPIB^N)
    ]

    #tpulse = induction_time .+ [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

    mdl = ODESystem(eqns, 
    continuous_events = vcat(
        [[t ~ induction_time] => [β1 ~ beta2, qSw ~ qSset2]],
        #[[t ~ induction_time] => [qSw ~ qSset2]],
        [[t ~ 7] => [qSw ~ qSset]]
        ); 
    name=name)

    return complete(structural_simplify(mdl))
end