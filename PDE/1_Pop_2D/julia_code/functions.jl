function dRdtFunc(R, V, kappaV, Delta, tau)
    return (1/tau)*(-kappaV * R .+ 2 * R .* V .+ Delta / (tau * pi))
end

function dVdtFunc(R, V, g, kappaS,tau, eta_0)
    return  (1/tau)*(kappaS*g .- (pi * tau * R) .^ 2 .+ V .^ 2 .+ eta_0)
end

function d2Rdt2Func(R, V, g, tau, Delta, kappaS, kappaV, eta_0)
    return (1 / tau) .* (-kappaV .+ 2 * V) .* F1(R, V, kappaV, Delta, tau) .+
           (2 / tau) .* R .* G1(R, V, g,kappaS, tau, eta_0)
end

function init_conds(nterms,kc,phi_vec,c_vec,X,Y)
    p = 0
    for i = 1:nterms
        k1,k2 = kc*[cos(phi_vec[i]) sin(phi_vec[i]) ]
        p = p .+ c_vec[i]*exp.(im*k1*X .+ im*k2*Y) + conj(c_vec[i])*exp.(-im*k1*X .- im*k2*Y)
    end
    return p
end
