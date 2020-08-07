using Plots
using SparseArrays
using ProgressMeter
using DifferentialEquations
using NPZ
using JLD
#include("psi_pde_neumannBC.jl")7
include("banded_matrices.jl")
include("RHS_funcs_wizardHat.jl")
include("functions.jl")
include("find_SS.jl")

#Editable parameters
#Spatial discretisation
X_max = 4pi #m
X_min = -4pi
N =  2^7
T_max = 50
#params
v = 2
eta_0=6
Delta = 0.5
alfa = 5.0
kappaV = 0.3
kappaS = 12.0
tau = 1
tau_a = 0.5
delta = 0.0


#Spatially discretisation
#
dx = 2*(X_max-X_min)/N

dxdx = dx * dx;
X=N
XX = Int(X^2)
Xspace =  LinRange(X_min,X_max,X)

zeroRV = SteadyState()
print("\n Fixed points are : R = ", zeroRV[1],"\n V = ", zeroRV[2],"\n")

R0 = zeros(XX)
V0=zeros(XX)
A10=zeros(XX)
A20=zeros(XX)
A30=zeros(XX)
A40=zeros(XX)
psi0 = zeros(XX)
g0 = zeros(XX)
p0=zeros(XX)
a0 = zeros(XX)
R0 = reshape(R0,X,X)

perturb = zeros(X,X)
perturb =  real(0.01*init_conds(1,3,[0,2pi/3,4pi/3],[1,1,1],Xspace,Xspace'))
sig = 2
guass_init = exp.(-(Xspace.^2)/sig .-((Xspace').^2)/sig)
perturb = reshape(guass_init,XX)

R0 = zeroRV[1] .+ perturb
V0 = zeroRV[2] .+ perturb
psi0 = perturb
A10 = perturb
A20 = perturb
A30 = perturb
g0 = perturb
p0 = perturb
a0 = perturb


#R0=5*rand(XX)
print("\n Use end of last simulation as initial conditions? (y/N)")
ans1 = readline(stdin)
if ans1 == "Y"
    #u0 = sol[:,end]
    u0 = zeros(9 * XX)
    #u0 = load("endsim.jld","arr")
    u0 = sol[:,end]
    Dxx = (1/dxdx)*D2r1(X, XX)
    Dyy = (1/dxdx)*D2r2(X, XX)
    Dxxyy = Dxx+Dyy
    #print(size(u0))
elseif ans1 == "y"
    u0 = sol[:,end]
else
print("\n Making sparse operator matrices ...... \n")
Dxx = (1/dxdx)*D2r1(X, XX)
Dyy = (1/dxdx)*D2r2(X, XX)
Dxxyy = Dxx+Dyy
print("done")
u0 = zeros(9 * XX)
u0[1:XX] .= R0
u0[XX + 1:2 * XX] .= V0
u0[2*XX + 1:3*XX] .= psi0
u0[3*XX + 1:4*XX] .= A10
u0[4*XX + 1:5*XX] .= A20
u0[5*XX + 1:6*XX] .= A30
u0[6*XX + 1:7*XX] .= p0
u0[7*XX + 1:8*XX] .= g0
u0[8*XX+1:9*XX] .= a0
end


#readline(stdin)u0 .= sol[:,end]
#u0 .= 5*rand(8*XX)





tspan = (0.0, T_max)
params = [tau,tau_a,kappaV,kappaS, eta_0, alfa, Delta, v , Dxxyy,delta]
prob = ODEProblem(rhsFun, u0, tspan, params)
print("\n Solving....... \n")
sol = solve(prob,saveat =.10, progress = true)
print("\n Done! \n")


R = reshape(sol[1:XX,:],X,X,length(sol))
V = reshape(sol[XX+1:2*XX,:],X,X,length(sol))
W = pi*R .+ im*V
Z = (1 .-conj.(W))./(1 .+conj.(W))
sync = abs.(Z)


print("\n Save time series? (y/N)")
ans = readline(stdin)
if ans == "Y" || ans == "y"
    print("Saving data...")
npzwrite("/home/james/PhD_Work/Python_Code/Brain_Top_Paper/Movie_Scripts/data/2D_planarwaves.npy", Z[:,:,:])
end


@gif for i = 1:1:size(sync,3)
heatmap(abs.(Z[:,:,i]))
end
