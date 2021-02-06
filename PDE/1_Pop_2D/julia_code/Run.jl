using Plots
using SparseArrays
using ProgressMeter
using DifferentialEquations
using NPZ
using JLD
# # # # # # # # # # # # # # # # # # # # # #
#
# Main file for running the model in
# 2 spatial dimension.
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
X_max = 1.0pi  # maximum of domain
X_min = -1.0pi # minimum of domain
X = 50       #number of grid points
T_max = 20000     # number of time points (dt is decided by solver)
dx = 2*(X_max-X_min)/X # spatial discretisation size
print("\n dx = ",dx,"\n")
dxdx = dx * dx;
XX = Int(X^2)
Xspace =  LinRange(X_min,X_max,X) # the spatial grid
# # # # # # # # # # # # # # # # # # # # # # # # # # #
#    MODEL PARAMETERS
#
v = 1 # axonal velocity
eta_0 =2 # mean drive
Delta = 0.5        # coherence
alfa = 0.5              # synaptic time constant
kappaV =0.695# gap junction strength
kappaS =12     # synaptic coupling strength
tau = 20# membrane time constant
saveat =50
# # # # # # # # # # # # # # # # # # # # # # # # # # #

# travelling & target waves kappaV =0.693
#
#
#

zeroRV = SteadyState()
print("\n Fixed points are : R = ", zeroRV[1],"\n V = ", zeroRV[2],"\n")
print("\n Initialising ... \n")
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

perturb = real(0.001*init_conds(3,1,[0,2pi/3,4pi/3],[1,1,1],Xspace,Xspace'))
perturb = real(0.001*init_conds(2,1,[0,2pi/3,4pi/3],[1,1,1],Xspace,Xspace'))
#perturb = real(0.001*init_conds(1,1,[0,2pi/3,4pi/3],[1,1,1],Xspace,Xspace')) #
#perturb = 0.001*exp.(-abs.(Xspace .+ Xspace'))
#perturb = 0.001*(exp.(-Xspace.^2 .- (Xspace').^2)) # guassian for target waves pattern
#perturb = real(0.001*init_conds(1,1,[0,2pi/3,4pi/3],[1,1,1],Xspace,Xspace')) #
#perturb = 0.01*cos.(1*Xspace)*cos.(1*Xspace')
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
        u0 = load("endsim2.jld","arr")
        #u0 = sol[:,end]
        Dxx = (1/dxdx)*D2r1(X, XX)
        Dyy = (1/dxdx)*D2r2(X, XX)
        Dxxyy = Dxx+Dyy
        #print(size(u0))
    else
        u0 = zeros(8 * XX)
        u0[1:XX] .= R0
        u0[XX + 1:2 * XX] .= V0
        u0[2*XX + 1:3*XX] .= psi0
        u0[3*XX + 1:4*XX] .= A10
        u0[4*XX + 1:5*XX] .= A20
        u0[5*XX + 1:6*XX] .= A30
        u0[6*XX + 1:7*XX] .= p0
        u0[7*XX + 1:8*XX] .= g0
        print("\n Making sparse operator matrices... ")
        Dxx = (1/dxdx)*D2r1(X, XX) # construct the differential operator matrices in x direction
        Dyy = (1/dxdx)*D2r2(X, XX) # construct the differential operator matrices in y direction
        Dxxyy = Dxx+Dyy # combine the differential operator matrices.
        print("DONE! \n")
    end




#solve
tspan = (0.0, T_max)
params = [tau,kappaV,kappaS, eta_0, alfa, Delta, v , Dxxyy]
prob = ODEProblem(rhsFun, u0, tspan, params)
print("\n Solving... \n")
sol = solve(prob,saveat =saveat, progress = true)
print("\n Done! \n")


R = reshape(sol[1:XX,:],X,X,length(sol))
Vc = reshape(sol[XX+1:2*XX,:],X,X,length(sol))
W = tau*pi*R .+ im*Vc
Z = (1 .-conj.(W))./(1 .+conj.(W))
t=LinRange(0,T_max,length(sol[1,:]))
sync = abs.(Z[:,:,end-100:end])


print("\n Save time series? (y/N)")
ans = readline(stdin)
if ans == "Y" || ans == "y"
    print("Saving data...")
    npzwrite("/home/james/PhD_Work/Python_Code/Brain_Top_Paper/Movie_Scripts/data/2D_standingwaves.npy", Z[:,:,:])
    print("Data saved.")
end


@gif for i = 1:1:size(sync,3)
heatmap(sync[:,:,i])
end

#save("/home/james/PhD_Work/Julia_Code/BrainTopPaper/PDE/1_Pop_2D/julia_code/endsim2.jld","arr", sol[:,end])
