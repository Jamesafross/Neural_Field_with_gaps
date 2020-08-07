using Plots
using SparseArrays
using ProgressMeter
using DifferentialEquations
using NPZ
include("banded_matrices.jl")
include("RHS_funcs_wizardHat.jl")
include("functions.jl")
include("find_SS.jl")

X_max = 2pi;
dx = pi/(2^8)
T_max = 1000
dxdx = dx * dx;
X = Int(2*X_max / dx);
#params
v = 2
eta_0 =-0.1
Delta = 0.5
alfa = 3
kappaV =1
kappaS = -15
tau = 1
k_c = 1
Dxx = (1 / dxdx) * D2x(X)
X_space = LinRange(-X_max,X_max,X)
print(X)
zeroRV = SteadyState()
print("\n Fixed points at: \n R = ", zeroRV[1], " \n V = ", zeroRV[2],"\n")

R0 = zeros(X)
V0=zeros(X)
A10=zeros(X)
A20=zeros(X)
A30=zeros(X)
A40=zeros(X)
psi0 = zeros(X)
g0 = zeros(X)
p0=zeros(X)
perturb =  .01*cos.(k_c * X_space)
R0 .= zeroRV[1] .+ 1*perturb
V0 .= zeroRV[2] .+ 1*perturb
psi0 .= 0 .+ 1*perturb
A10 .= 0 .+ 1*perturb
A20 .= 0 .+ 1*perturb
p0 .= 0 .+ 1*perturb
g0 .= 0 .+ 1*perturb
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

W = pi*R + im*V
Z = (1 .-conj.(W)) ./(1 .+conj.(W))
sync = abs.(Z)
syncSave = sync[:,end-200:end]
RSave = R[:,end-200:end]

#@gif for i = 1:length(syncSave[1,:])
#  plot(syncSave[:,i],ylims=[minimum(syncSave),maximum(syncSave)])
#end
print(size(syncSave))

print("save? (Y/N) ")
ans = readline()
if ans == "Y"
npzwrite("/home/james/PhD_Work/Python_Code/Brain_Top_Paper/Review_Paper/data/1D_sync_Turing.npy", syncSave)
end

heatmap(syncSave)
@gif for i = 1:2000; plot(sync[:,i],ylims=[0,1]);end