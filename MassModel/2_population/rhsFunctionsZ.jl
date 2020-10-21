function V(z)
    return imag((1 - conj.(z)) / (1 + conj.(z)))
end

function FR(z,tau)
    return (tau/pi)*((1-abs(z).^2)./(1+z+conj(z)+abs(z).^2))
end

function F(z, eta0, delta, kappa)
    return 0.5 * (-im * (z - 1)^2 +
           (z + 1)^2 * (-delta + im * eta0))

end

function G(z, g1, g2, v_syn1, v_syn2)
    return 0.5 * (im * ((z + 1)^2) * g1 * v_syn1 - (z^2 - 1) * g1) +
           0.5 * (im * ((z + 1)^2) * g2 * v_syn2 - (z^2 - 1) * g2)

end



function rhs_functionZ(du, u, p, t)
    #zE
    du[1] =
        (1/tauE)*(F(u[1], eta_0E, deltaE, kappaVE) + G(u[1], u[4], u[8], VsynEE, VsynEI))
    #zI
    du[2] =
        (1/tauE)*(F(u[2], eta_0I, deltaI, kappaVI) + G(u[2], u[10], u[6], VsynII, VsynIE))
    #psiEE
    du[3] = alphaEE * (-u[3] + betaEE * FR(u[1],tauE))
    #gEE
    du[4] = alphaEE * (-u[4] + u[3])
    #psiIE
    du[5] = alphaIE * (-u[5] + betaIE * FR(u[1],tauE))
    #gIE
    du[6] = alphaIE * (-u[6] + u[5])
    #psiEI
    du[7] = alphaEI * (-u[7] + betaEI * FR(u[2],tauI))
    #gEI
    du[8] = alphaEI * (-u[8] + u[7])
    #psiII
    du[9] = alphaII * (-u[9] + betaII * FR(u[2],tauI))
    #gII
    du[10] = alphaII * (-u[10] + u[9])

end
