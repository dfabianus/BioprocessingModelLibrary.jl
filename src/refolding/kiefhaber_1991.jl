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

Kiefhaber_network = @reaction_network begin
    k_n, I --> N
    k_a, 2*I --> A
  end