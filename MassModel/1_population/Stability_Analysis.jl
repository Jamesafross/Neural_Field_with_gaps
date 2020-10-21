using Roots
using LinearAlgebra
include("findSteadyState.jl")
include("jacobian.jl")

v_syn = -20
eta_0 =1;
delta = 0.2;
Alpha = 0.4;
beta = 1;
a=0.01
b=1

N=1000
evals = complex(zeros(4,N))
kappa_vec = LinRange(0,2,N)
for i = 1:N
global kappa = kappa_vec[i]
p = [v_syn,kappa,eta_0,delta,Alpha]
ss = (SS(p,a,b))
u1 = ss[1]
u2 = ss[2]
u3 = ss[3]
evals[:,i] = eigen(J(u1,u2,u3,kappa,v_syn,beta,Alpha)).values
end
#print(evals)

#scatter(evals)
