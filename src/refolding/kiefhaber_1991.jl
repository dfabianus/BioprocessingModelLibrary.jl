k_n_fun(a_n, b_n, D) = a_n * (1 + D) ^ b_n
k_a_fun(a_a, b_a, D) = a_a * (1 + D) ^ b_a

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

@mtkmodel Kiefhaber_B begin
    @parameters begin
        n = 2 # Order of aggregation
        a_n = 1.3343
        b_n = -8.6824
        a_a = 12.0465
        b_a = -16.7869
        D = 0.1
    end
    @variables begin
        I(t)
        N(t)
        A(t)
        k_n(t) # refolding rate
        k_a(t) # aggregation rate
    end
    @equations begin
        Dt(I) ~ -k_n*I - k_a*I ^ n 
        Dt(N) ~ k_n*I
        Dt(A) ~ k_a*I ^ n
        k_n ~ k_n_fun(a_n, b_n, D)
        k_a ~ k_a_fun(a_a, b_a, D)
    end
end

function Kiefhaber_network()
    Kiefhaber_network = @reaction_network begin
        k_n_fun(a_n, b_n, D), I --> N
        k_a_fun(a_a, b_a, D), 2*I --> A
    end
    return Kiefhaber_network
end

# requires Graphviz
#Graph(Kiefhaber_network())