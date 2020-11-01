using Plots
using SparseArrays
using ProgressMeter
using DifferentialEquations
using NPZ
using JLD
# # # # # # # # # # # # # # # # # # # # # #
#
# Main file for running the model in
# 1 spatial dimension.
#
# You are able to edit parameters directly
# in this file.
#
# # # # # # # # # # # # # # # # # # # # # #

# include files
include("banded_matrices.jl")
include("RHS_funcs_wizardHat.jl")
include("functions.jl")
include("find_SS.jl")

# # # # # # # # # # # # # # # # # # # # # # # # # # #
# SPATIAL DISCRETISATION PARAMETERS
#
X_max = pi # maximum of domain
X_min = -pi # minimum of domain
X = 200 # number of grid points
T_max =100 # number of time points (dt is decided by adaptive solver)
dx = 2*(X_max-X_min)/X # spatial discretisation size
dxdx = dx * dx;
XX = Int(X^2)
Xspace =  LinRange(X_min,X_max,X) # the spatial grid
Dxx = (1/dxdx)*D2r1(X, XX) # construct the differential operator matrices in x direction
Dyy = (1/dxdx)*D2r2(X, XX) # construct the differential operator matrices in y direction
Dxxyy = Dxx+Dyy # combine the differential operator matrices.
# # # # # # # # # # # # # # # # # # # # # # # # # # #
#    MODEL PARAMETERS
#
v = 3                 # axonal velocity
eta_0 =0.3            # mean drive
Delta = 0.5           # coherence
alfa = 3              # synaptic time constant
kappaV =1.2           # gap junction strength
kappaS = 5            # synaptic coupling strength
tau = 1               # membrane time constant
saveat = .1
# # # # # # # # # # # # # # # # # # # # # # # # # # #


zeroRV = SteadyState()
print("\n Fixed points are : R = ", zeroRV[1],"\n V = ", zeroRV[2],"\n")

# Initialising initial conditions
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
gauss_init = exp.(-Xspace.^2 .- (Xspace').^2)
perturb =  real(0.01*init_conds(3,3,[pi/2,2pi/3,4pi/3],[1,3,3],Xspace,Xspace'))
#perturb = gauss_init
#perturb = 0.01*cos.(3*Xspace)*cos.(3*Xspace')

perturb = reshape(perturb,XX)

# initial conditions
R0 = zeroRV[1] .+ perturb
V0 = zeroRV[2] .+ perturb
psi0 = perturb
A10 = perturb
A20 = perturb
A30 = perturb
g0 = perturb
p0 = perturb

print("\n Use end of last simulation as initial conditions? (y/N)")
ans1 = readline(stdin)
if ans1 == "Y" || ans1 == "y"
    #u0 = sol[:,end]
    u0 = zeros(9 * XX)
    #u0 = load("endsim.jld","arr")
    u0 = sol[:,end]
    Dxx = (1/dxdx)*D2r1(X, XX)
    Dyy = (1/dxdx)*D2r2(X, XX)
    Dxxyy = Dxx+Dyy
    #print(size(u0))
else
print("\n Making sparse operator matrices ...... \n")
print("done")
u0 = zeros(8 * XX)
u0[1:XX] .= R0
u0[XX + 1:2 * XX] .= V0
u0[2*XX + 1:3*XX] .= psi0
u0[3*XX + 1:4*XX] .= A10
u0[4*XX + 1:5*XX] .= A20
u0[5*XX + 1:6*XX] .= A30
u0[6*XX + 1:7*XX] .= p0
u0[7*XX + 1:8*XX] .= g0
end

#solve
tspan = (0.0, T_max)
params = [tau,kappaV,kappaS, eta_0, alfa, Delta, v , Dxxyy]
prob = ODEProblem(rhsFun, u0, tspan, params)
print("\n Solving....... \n")
sol = solve(prob,saveat =saveat, progress = true)
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
npzwrite("/home/james/PhD_Work/Python_Code/Brain_Top_Paper/Movie_Scripts/data/2D_TuringHopf_kappaV=3.npy", Z[:,:,:])
end


@gif for i = 1:1:size(sync,3)
heatmap(abs.(Z[:,:,i]))
end
