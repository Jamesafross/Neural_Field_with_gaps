function F1(R, V, kappaV, Delta, tau)
    return (1/tau)*(-kappaV * R .+ 2 * R .* V .+ Delta / (tau * pi))
end

function G1(R, V, g, tau, eta_0)
    return  (1/tau)*(g .- (pi * tau * R) .^ 2 .+ V .^ 2 .+ eta_0)
end

function ddR(R, V, g, tau, Delta, kappaV, eta_0)
    return (1 / tau) .* (-kappaV .+ 2 * V) .* F1(R, V, kappaV, Delta, tau) .+
           (2 / tau) .* R .* G1(R, V, g,tau, eta_0)
end


function dddR(R, V, g, kappa_v, tau, Delta, eta_0)
    return (1 / tau^2) .* (
        (.-kappa_v .+ 2 .* V) .^ 2 .+
        (2 / tau^2) .*
        ( g .- 3 * (pi * tau * R) .^ 2 .+ V .^ 2 .+ eta_0)
    ) .* F1(R, V, kappa_v, Delta, tau) .+
           (1 / tau^2) .* (
        .-2 .* kappa_v .* R .+
        2 .* (-kappa_v .* R .+ 4 .* R .* V .+ Delta / (tau * pi)) .+
        4 .* R .* V
    ) .* G1(R, V, g, tau, eta_0)
end
