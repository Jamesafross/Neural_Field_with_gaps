using Plots
using SparseArrays
using ProgressMeter
using DifferentialEquations
using NPZ
#Include Files
include("banded_matrices.jl")
include("RHS_funcs_wizardHat.jl")
include("functions.jl")
include("find_SS.jl")
#spatial parameters
X_max = 8pi;
dx = pi/(2^6)
T_max = 600
dxdx = dx * dx;
X = Int(2*X_max / dx);
Dxx = (1 / dxdx) * D2x(X)
X_space = LinRange(-X_max,X_max,X)
#model parameters
v = 3    # axonal velocity
eta_0 =0.3 # mean drive
Delta = 0.5 # coherence
alfa = 3 # synaptic time constant
kappaV =1.2 # gap junction strength
kappaS = 5 #synaptic coupling strength
tau = 1 # membrane time constant


k_c = 1 #wave length to excite (initial conditions)

print(X)
#steady state
zeroRV = SteadyState()
print("\n Fixed points at: \n R = ", zeroRV[1], " \n V = ", zeroRV[2],"\n")

#initialise initial conditions
R0 = zeros(X)
V0=zeros(X)
A10=zeros(X)
A20=zeros(X)
A30=zeros(X)
A40=zeros(X)
psi0 = zeros(X)
g0 = zeros(X)
p0=zeros(X)
#perturbation to the steady state
perturb =  0.0001*cos.(X_space)
#initial conditions
R0 .= zeroRV[1] .+ 1*perturb
V0 .= zeroRV[2] .+ 1*perturb
psi0 .= 0 .+ 1*perturb
A10 .= 0 .+ 1*perturb
A20 .= 0 .+ 1*perturb
p0 .= 0 .+ 1*perturb
g0 .= 0 .+ 1*perturb

#initial conditions of the solver
u0 = zeros(8 * X)
u0[1:X] .= R0
u0[X + 1:2 * X] .= V0
u0[2*X + 1:3*X] .= psi0
u0[3*X + 1:4*X] .= A10
u0[4*X + 1:5*X] .= A20
u0[5*X + 1:6*X] .= A30
u0[6*X + 1:7*X] .= p0
u0[7*X + 1:8*X] .= g0

#solve
tspan = (0.0, T_max)
params = [tau, kappaV,kappaS, eta_0, alfa, Delta, v , Dxx]
prob = ODEProblem(rhsFun, u0, tspan, params)
print("Solving...")
sol = solve(prob,saveat = 0.2, progress = true)
print("\n Done!")

R = sol[1:X,:]
V = sol[X + 1:2*X,:]
g = sol[7*X+1:8*X,:]

#conversion to kuromoto frame
W = pi*R + im*V
Z = (1 .-conj.(W)) ./(1 .+conj.(W))
sync = abs.(Z)
syncSave = sync[:,end-200:end]
RSave = R[:,end-200:end]


#save file, plotting, etc
#@gif for i = 1:length(syncSave[1,:])
#  plot(syncSave[:,i],ylims=[minimum(syncSave),maximum(syncSave)])
#end
print(size(syncSave))

#print("save? (Y/N) ")
#ans = readline()
#if ans == "Y"
#npzwrite("/home/james/PhD_Work/Python_Code/Brain_Top_Paper/Review_Paper/data/1D_sync_Turing.npy", syncSave)
#end

heatmap(syncSave[:,1:end])
#@gif for i = 1:2000; plot(sync[:,i],ylims=[0,1]);end
