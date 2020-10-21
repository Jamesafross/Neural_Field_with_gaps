using FFTW
using Plots
using Statistics
using LinearAlgebra
using NPZ
using DSP
using MAT
using SparseArrays
using DifferentialEquations
include("Solvers.jl")
include("trim.jl")
include("rhsFunctions.jl")
include("rhsFunctionsZ.jl")
include("firingrate.jl")

#parameters
kappaV = 1.2;
kappaS=1.0;
tau=1
eta_0 = 1;
delta = 0.5;
Alpha = 1;
beta = 1;




#time discretisation
T = 100

dt = 0.1


u01 = zeros(4)
u02 = complex(zeros(3))
u01 = rand(4)
saveat = dt
p = [kappaV,kappaS,tau,eta_0,delta,Alpha]

tspan = (0.0,1200)
prob1 = ODEProblem(rhs_function,u01,tspan,p)
sol1 = solve(prob1,reltol=1e-8, abstol=1e-8)





Fs = 20000




npzwrite("/home/james/PhD_Work/Python_Code/Brain_Top_Paper/Review_Paper/data/massmodel_periodic_onepop.npy", sol1[:,:])
plot(sol1[1,end-1000:end])
