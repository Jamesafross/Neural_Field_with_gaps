using Plots
using SparseArrays
using ProgressMeter
using DifferentialEquations
using NPZ
# # # # # # # # # # # # # # # # # # # # # #
#
# Main file for running the model in
# 1 spatial dimension.
#
# You are able to edit parameters directly
# in this file.
#
# # # # # # # # # # # # # # # # # # # # # #

#Include Files
include("banded_matrices.jl")
include("RHS_funcs_wizardHat.jl")
include("functions.jl")
include("find_SS.jl")
# # # # # # # # # # # # # # # # # # # # # # # # # # #
# SPATIAL DISCRETISATION PARAMETERS

X_max = 5pi; # size f domain
dx = pi/(2^5) # spatial discretisation size
T_max = 100000 # maximum time
dxdx = dx * dx;
X = Int(2*X_max / dx); # number of grid points
Dxx = (1 / dxdx) * D2x(X) # create second order finite difference matrix
X_space = LinRange(-X_max,X_max,X)

# # #    MODELPARAMETERS  # # # # # # # # # # # # # #

v = .11       # axonal velocity
eta_0 =1   # mean drive
Delta = 0.5   # coherence
alfa = 0.5     # synaptic time constant
kappaV =  1.8# gap junction strength
kappaS =10    # synaptic coupling strength
tau = 15      # membrane time constant
# # # # # # # # # # # # # # # # # # # # # # # # # # #


#TH1: v = 0.11, kappaV = 0.855
#TH2: v = 1, kappaV = 0.88
#H: v=0.1 kappaV = 0.85

#from paper; v = 0.11
#low kappaV = 0.855
#mid kappaV = 0.864
#high kappaV = 1.0
#very high kappaV = 1.2

dt=0.5
k_c = 1 # wave length to excite (initial conditions)

print(X)
#steady state
zeroRV = SteadyState()
W_SS = pi*tau*zeroRV[1] + im*zeroRV[2]
sync_SS = abs.((1 .-conj.(W_SS)) ./(1 .+conj.(W_SS)))
print("\n Fixed points at: \n R = ", zeroRV[1], " \n V = ", zeroRV[2],"\n"," \n |Z| = ", sync_SS,"\n")

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
perturb =  0.001*cos.(k_c*X_space)
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
sol = solve(prob,saveat = dt, progress = true,maxiters=1e7)
print("\n Done!")

R = sol[1:X,:]
V = sol[X + 1:2*X,:]
g = sol[7*X+1:8*X,:]

#conversion to kuromoto frame
W = pi*tau*R + im*V
Z = (1 .-conj.(W)) ./(1 .+conj.(W))
sync = abs.(Z)
syncSave = sync[:,end-500:end]
RSave = R[:,end-500:end]
t = collect(0:dt:T_max




#save file, plotting, etc
#@gif for i = 1:length(syncSave[1,:])
#  plot(syncSave[:,i],ylims=[minimum(syncSave),maximum(syncSave)])
#end

print("save? (Y/N) ")
ans = readline()
if ans == "Y" || ans == "y"
npzwrite("/home/james/PhD_Work/Python_Code/Brain_Top_Paper/Review_Paper/data/1D_sync_highhighk.npy", syncSave)
end

M = 500
heatmap(t[end-M:end] .-t[end-M],X_space,sync[:,end-M:end])
#@gif for i = 1:2000; plot(sync[:,i],ylims=[0,1]);end
