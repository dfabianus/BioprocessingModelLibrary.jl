@mtkmodel Kiefhaber begin
    @parameters begin
        n = 2
        k_n = 0.1
        k_a = 0.1
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