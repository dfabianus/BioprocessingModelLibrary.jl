using ModelingToolkit

Am = [1 0 3; 2 0 1; 0 1 4]
Bm = [4 5; 3 2; 1 4]

@parameters t
Dt = Differential(t)

@assert size(Am, 1) == size(Bm, 1)
nr_of_components = size(Am, 1) # 3
nr_of_incoming_feeds = size(Am, 2) # 3
nr_of_independent_rates = size(Bm, 2) # 2

@parameters begin
    A[1:nr_of_components, 1:nr_of_incoming_feeds]
    B[1:nr_of_components, 1:nr_of_independent_rates]
end

@variables begin
    u_1(t)[1:nr_of_incoming_feeds] # [2 x 1] vector
    c(t)[1:nr_of_components] # [3 x 1] vector
    q_i(t)[1:nr_of_independent_rates] # [2 x 1] vector
end

eqns = [
    Dt.(c) ~ A * u_1 + B * q_i * c[1] 
]
# This gives:
# 1-element Vector{Symbolics.Arr{Any, 1}}:
# (broadcast(~, broadcast(Differential(t), c(t)), broadcast(+, A*u_1(t), B*broadcast(*, q_i(t), Ref((c(t))[1])))))[Base.OneTo(3)]

@named mdl = ODESystem(eqns)
