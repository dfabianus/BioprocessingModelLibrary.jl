k_m_fun(a_m, b_m, D) = a_m * (1 + D) ^ b_m

@mtkmodel INAM begin
    @extend Kiefhaber_B()
    @parameters begin
        a_m
        b_m
    end
    @variables begin
        M(t)
        k_m(t)
    end
    @equations begin
        Dt(M) ~ k_m*I
        k_m ~ k_m_fun(a_m, b_m, D)
    end
end

function INAM_network()
    INAM_network = @reaction_network begin
        k_n_fun(a_n, b_n, D), I --> N
        k_a_fun(a_a, b_a, D), 2*I --> A
        k_m_fun(a_m, b_m, D), I --> M
    end
    return INAM_network
end