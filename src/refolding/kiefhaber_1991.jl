@mtkmodel Kiefhaber begin
    @parameters begin
        n = 2 # Order of aggregation
        k_n = 0.1 # refolding rate
        k_a = 0.1 # aggregation rate
    end

    @variables begin
        I(t)
        N(t)
        A(t)
    end

    @equations begin
        Dt(I) ~ -k_n*I - k_a*I ^ n 
        Dt(N) ~ k_n*I
        Dt(A) ~ k_a*I ^ n
    end
end

function Kiefhaber_network()
    k_n(a_n, b_n, D) = a_n * (1 + D) ^ b_n
    k_a(a_a, b_a, D) = a_a * (1 + D) ^ b_a

    Kiefhaber_network = @reaction_network begin
        k_n(a_n, b_n, D), I --> N
        k_a(a_a, b_a, D), 2*I --> A
    end

    return Kiefhaber_network
end