module OControllerProx;
using JuMP;
using MathProgBase;
using ProxSDP;
using LinearAlgebra;

export solveControllerObserver;

function solveControllerObserver(Ax,Bx,Cx,gamma,Delta1,Delta2,nx,mx,px)

#ky = Model(solver = ProxSDP.Optimizer());
# CONTROLLER DESIGN
ky=Model(with_optimizer(ProxSDP.Optimizer, log_verbose=false)); # to find control gain
#set_silent(ky)
#MOI.set(ky, MOI.Silent(), true);

    @variable(ky, X[1:nx,1:nx], PSD) ;
    @variable(ky, Y[1:nx,1:mx]) ;
    @variable(ky,o_k);

@constraint(ky, o_k==1);
@SDconstraint(ky, X >= 0.001*Matrix{Float64}(I, nx, nx));
@SDconstraint(ky, X*Ax'+Y*Bx'+Ax*X+Bx*Y'+gamma*X <= 0);

@objective(ky, Min, o_k); # feasibility problem

JuMP.optimize!(ky);

#status = solve(ky);

JuMP.objective_value(ky);

Xv=JuMP.value.(X);
Yv=JuMP.value.(Y);

P=pinv(Xv);
K=Yv'*P';


# OBSERVER DESIGN

ly=Model(with_optimizer(ProxSDP.Optimizer, log_verbose=false)); # to find observer gain

#MOI.set(ly, MOI.Silent(), true)
@variable(ly, Q[1:nx,1:nx], PSD) ;
@variable(ly, Z[1:px,1:nx]) ;

@variable(ly,o_l);

@constraint(ly, o_l==1);
@SDconstraint(ly, Q >= 0.001*Matrix{Float64}(I, nx, nx));
@SDconstraint(ly, Ax'*Q-Cx'*Z+Q*Ax-Z'*Cx+2*gamma*Q <= 0);

@objective(ly, Min, o_l); # feasibility problem

JuMP.optimize!(ly);

#status = solve(ky);

JuMP.objective_value(ly);

Qv=JuMP.value.(Q);
Zv=JuMP.value.(Z);

Q=Qv;
L=pinv(Qv)'*Zv';

# OUTPUT EXTRACT

return(P,K,Q,L);


end

end
