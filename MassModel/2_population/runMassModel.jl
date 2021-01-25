using FFTW
using Plots
using Statistics
using LinearAlgebra
using NPZ
using DSP
using MAT
using SparseArrays
using DifferentialEquations
include("trim.jl")
include("rhsFunctions.jl")
include("ICs.jl")

#parameters
tauE = 1
tauI= 1
kappaVE = 0.5
kappaVI = 0.5
kappaSEE = 15
kappaSIE = 25
kappaSEI = -15
kappaSII = -15
eta_0E = 5
eta_0I = -3
deltaE = 0.5
deltaI = 0.5
alphaEE = 0.2
alphaIE = 0.1
alphaEI = 0.07
alphaII = 0.06




#time discretisation
T = 10000
saveat = .1


# Initial Conditions
# from the ICs.jl file
u01= hopf_init()


p = [tauE,tauI,kappaVE,kappaVI,kappaSEE,kappaSIE,kappaSEI,kappaSII,eta_0E,eta_0I,deltaE,deltaI,alphaEE,alphaIE,alphaEI,alphaII]

tspan = (0.0,T)
prob1 = ODEProblem(rhs_function,u01,tspan,p)
sol1 = solve(prob1,saveat = 0.1)

#plot(real.(Z_to_W.(sol3[1,1:end]))/pi)
t = LinRange(0,T,length(sol1[1,:]))
npzwrite("/home/james/PhD_Work/Python_Code/Brain_Top_Paper/Review_Paper/data/massmodel_periodic_twopop_kappa=0.5.npy", sol1[:,end-5000:end])
plot(t[end-5000:end].-t[end-5000],[sol1[1,end-5000:end],sol1[2,end-5000:end]])
