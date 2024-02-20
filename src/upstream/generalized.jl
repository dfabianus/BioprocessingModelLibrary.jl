### d(c*V)/dt =  A * q_i * c_X * V + B * u1 - sum(u2) * c
using ModelingToolkit

Am = [4 5; 3 2; 1 4]
Bm = [1 0 3; 2 0 1; 0 1 4]
#name = "BioreactorDynamics"

@parameters t
Dt = Differential(t)

#function BioreactorDynamics(Am, Bm; name)

    @assert size(Am, 1) == size(Bm, 1)
    nr_of_independent_rates = size(Am, 2) # 2
    nr_of_components = size(Am, 1) # 3
    nr_of_incoming_feeds = size(Bm, 2) # 3

    @parameters begin
        A[1:nr_of_components, 1:nr_of_independent_rates]
        B[1:nr_of_components, 1:nr_of_incoming_feeds]
    end

    @variables begin
        V(t) # scalar
        u_1(t)[1:nr_of_incoming_feeds] # [2 x 1] vector
        #u_2(t)[1:2] # [2 x 1] vector
        c(t)[1:nr_of_components] # [3 x 1] vector
        q_i(t)[1:nr_of_independent_rates] # [2 x 1] vector
    end

    eqns = [
        Dt.(V) ~ sum(u_1)
        Dt.(c) ~ B * u_1 / V + A * q_i * c[1] #- sum(u_2) * c / V
    ]

    #tpulse = induction_time .+ [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

    @named mdl = ODESystem(eqns, t; name=name)

    #return complete(structural_simplify(mdl))
#end