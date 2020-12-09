# functions used in the solver
# R1: ∂R/∂t
# V1: ∂V/∂t
# R2: ∂²R/∂t²
# R3: ∂³R/∂t³

function R1(R, V, kappaV, Delta, tau)
    return (1/tau)*(-kappaV * R .+ 2 * R .* V .+ Delta / (tau * pi))
end

function V1(R, V, g,kappaS, tau, eta_0)
    return  (1/tau)*(kappaS*g .- (pi * tau * R) .^ 2 .+ V .^ 2 .+ eta_0)
end

function R2(R, V, g,kappaS,kappaV,tau, Delta, eta_0)
    return (1 / tau) .* (-kappaV .+ 2 * V) .* R1(R, V, kappaV, Delta, tau) .+
           (2 / tau) .* R .* V1(R, V, g,kappaS,tau, eta_0)
end


function R3(R, V, g, kappaS, kappaV, tau, Delta, eta_0,g_dot)

    #H = ∂²R/∂t²
    #∂H/∂t = ∂H/∂R * ∂R/∂t +  ∂H/∂V * ∂V/∂t  + ∂H/∂g * ∂g/∂t

    G = V1(R, V, g, kappaS,tau, eta_0)
    F = R1(R, V, kappaV, Delta, tau)
    dHdR = (1/tau^2).*(2 .*V .- kappaV).^2 .+ (2/tau).*G .+
           (1/tau^2).*(2 .*R).*(-2 .*R.*pi.^2 .*tau.^2)

    dHdV = (2/tau).*F .+ (1/tau^2).*(2 .*V .- kappaV).*(2 .*R) .+
           (1/tau^2).*(2 .*R).*(2 .*V)

    dHdg = (1/tau^2).*(2 .*R).*(kappaS)

    return dHdR .*F .+ dHdV.*G .+ dHdg.*g_dot
end
