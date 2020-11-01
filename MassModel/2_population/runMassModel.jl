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
VsynEE=14
VsynIE=8
VsynEI=-8
VsynII=-12
tauE = 1
tauI= 0.5
kappaVE = 0.5
kappaVI = 0.5
kappaSEE = 2.5
kappaSIE = 14
kappaSEI = -5.1
kappaSII = -1.4
eta_0E =5
eta_0I = -6
deltaE = 0.5
deltaI = 0.5
alphaEE = 1
alphaIE = 1.4
alphaEI = 0.7
alphaII = 0.4
betaEE = 2.0
betaIE = 1.0
betaEI = 2.0
betaII = 3



#time discretisation
T = 1000
saveat = .1


u01= hopf_init()


p = [tauE,tauI,kappaVE,kappaVI,kappaSEE,kappaSIE,kappaSEI,kappaSII,eta_0E,eta_0I,deltaE,deltaI,alphaEE,alphaIE,alphaEI,alphaII,betaEE,betaIE,betaEI,betaII]
p2 = [tauE,tauI,kappaVE,kappaVI,VsynEE,VsynIE,VsynEI,VsynII,eta_0E,eta_0I,deltaE,deltaI,alphaEE,alphaIE,alphaEI,alphaII,betaEE,betaIE,betaEI,betaII]

tspan = (0.0,T)
prob1 = ODEProblem(rhs_function,u01,tspan,p)
sol1 = solve(prob1,RK4())

#plot(real.(Z_to_W.(sol3[1,1:end]))/pi)

#npzwrite("/home/james/PhD_Work/Python_Code/Brain_Top_Paper/Review_Paper/data/massmodel_periodic_twopop_kappa=0.5.npy", sol1[:,:])
plot([sol1[1,500:1000],sol1[2,500:1000]])
