using Plots
using NPZ
using SparseArrays
using DifferentialEquations
include("Solvers.jl")
include("trim.jl")
include("rhsFunctions.jl")
include("firingrate.jl")

#parameters
kappaV = 1.2;
kappaS=1.0;
tau=15
eta_0 = 1
delta = 0.5;
Alpha = 0.5;
beta = 1;




#time discretisation
T = 2000
dt = 0.01


u01 = zeros(4)
u02 = complex(zeros(3))
u01 = rand(4)
saveat = dt
p = [kappaV,kappaS,tau,eta_0,delta,Alpha]

tspan = (0.0,T)
prob1 = ODEProblem(rhs_function,u01,tspan,p)
sol1 = solve(prob1,saveat=saveat,reltol=1e-8, abstol=1e-8)

R = sol1[1,:]
V = sol1[2,:]
W = pi*R + im*V
Z = (1 .-conj.(W)) ./(1 .+conj.(W))
sync = abs.(Z)

tvec = collect(0:0.01:1000)
sol1_1ms = sol1[:,end-100000:end]
tvec_1ms = tvec[end-100000:end].-tvec[end-100000]
npzwrite("/home/james/PhD_Work/Python_Code/Brain_Top_Paper/Review_Paper/data/massmodel_periodic_onepop.npy", sol1_1ms)
plot(tvec_1ms, [sol1_1ms[2,:]])
