function pulseRefoldingModel(pulse_times, mP_pulses, mD_pulses, V_pulses; name)

    eqs = [
        Dt(V) ~ 0
        Dt(Ii) ~ -kₙ*Ii - kₐ*Ii ^ n 
        Dt(N) ~ kₙ*Ii
        Dt(A) ~ kₐ*Ii ^ n
        Dt(D) ~ 0
        cI  ~ Ii / V
        cN  ~ N / V
        cA  ~ A / V
        cD  ~ D / V
        cP  ~ P / V
        Y ~ Dt(N) / (Dt(N) + Dt(A))
        Y_cum ~ N / (N + A)
        kₙ ~ aₙ * (1+cD) ^ bₙ
        kₐ ~ aₐ * (1+cD) ^ bₐ
        P ~ Ii + N
        Intensity_mon ~ monod(cP, [kInt3, kInt4])
        delta_aew ~ -model_aew((cN+cA)/cP, [kAew1, kAew2, kAew3])
    ]

    times = pulse_times[2:end]
    mds = D .~ D .+ mD_pulses[2:end]
    mps = Ii .~ Ii .+ mP_pulses[2:end]
    Vs = V .~ V .+ V_pulses[2:end]

    sys = ODESystem(eqs, t,
    discrete_events = vcat(
        [[a] => [b] for (a,b) in zip(times, mds)],
        [[a] => [b] for (a,b) in zip(times, mps)],
        [[a] => [b] for (a,b) in zip(times, Vs)],
    ); name = name
    )
    return sys
end