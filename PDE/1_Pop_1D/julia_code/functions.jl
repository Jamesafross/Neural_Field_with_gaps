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

function R2(R, V, g, tau, Delta, kappaS, kappaV, eta_0)
    return (1 / tau) .* (-kappaV .+ 2 * V) .* R1(R, V, kappaV, Delta, tau) .+
           (2 / tau) .* R .* V1(R, V, g,kappaS, tau, eta_0)
end


function R3(R, V, g, kappaS, kappaV, tau, Delta, eta_0, g_dot)

    #H(R,V,g) = ∂²R/∂t²
    #∂³R/∂t³ = ∂H/∂t = ∂H/∂R * ∂R/∂t +  ∂H/∂V * ∂V/∂t  + ∂H/∂g * ∂g/∂t

    F = R1(R, V, kappaV, Delta, tau)
    G = V1(R, V, g, kappaS,tau, eta_0)

    dHdR = (1/tau^2).*(2 .*V .- kappaV).^2 .+ (2/tau).*G .+
           (1/tau^2).*(2 .*R).*(-2 .*R.*(pi^2) .*(tau^2))

    dHdV = (2/tau).*F .+ (1/tau^2).*(2 .*V .- kappaV).*(2 .*R) .+
           (1/tau^2).*(2 .*R).*(2 .*V)

    dHdg = (1/tau^2).*(2 .*R).*(kappaS)

    return dHdR .*F .+ dHdV.*G .+ dHdg.*g_dot


    #return (1/(pi*tau^4)).*(2 .*(Delta .- 2*pi*tau.*R.*(kappaV .- 3 .*V)).*
    #       (eta_0 .+ kappaS.*g - pi.^2*tau.^2 .*R.^2 .+ V.^2) .+ (Delta .- pi.*tau.*R.*(kappaV .- 2 .*V)).*
    #       (2 .*eta_0 .+ kappaV^2 .+ 2 .*kappaS.*g .- 6 .*pi^2*tau^2 .*R.^2 .- 4*kappaV .*V .+ 6 .*V.^2) +
    #       2*pi*kappaS*tau^2 .*R .*g_dot)
end



function dddR(R, V, g, kappaS, kappaV, tau, Delta, eta_0, g_dot)
    return (1 / tau^2) .* (
    (.-kappaV .+ 2 .* V) .^ 2 .+
        (2 / tau^2) .*
        ( g .- 3 * (pi * tau * R) .^ 2 .+ V .^ 2 .+ eta_0)
    ) .* R1(R, V, kappaV, Delta, tau) .+
           (1 / tau^2) .* (
        .-2 .* kappaV .* R .+
        2 .* (-kappaV .* R .+ 4 .* R .* V .+ Delta / (tau * pi)) .+
        4 .* R .* V
    ) .* V1(R, V, g, kappaS, tau, eta_0)
end
