function rhs_function(du, u, p, t)
    rE=u[1]
    rI=u[2]
    vE=u[3]
    vI=u[4]
    pEE=u[5]
    gEE=u[6]
    pIE=u[7]
    gIE=u[8]
    pEI=u[9]
    gEI=u[10]
    pII=u[11]
    gII=u[12]
    #rE
    du[1] =(1/tauE)*(-kappaVE * rE + 2 * rE * vE + (deltaE / (tauE*pi)))
    #rI
    du[2] =(1/tauI)*(-kappaVI * rI + 2 * rI * vI + (deltaI / (tauI*pi)))
    #vE
    du[3] =(1/tauE)*(gEE * kappaSEE + gEI * kappaSEI - (tauE^2)*(pi^2)*(rE^2) +  vE^2 + eta_0E)
    #vI
    du[4] =(1/tauI)*(gIE * kappaSIE + gII *  kappaSII - (tauI^2)*(pi^2)*(rI^2) + vI^2 + eta_0I)
    #psiEE
    du[5] = alphaEE * (-pEE +  rE)
    #gEE
    du[6] = alphaEE * (-gEE + pEE)
    #psiIE
    du[7] = alphaIE * (-pIE +  rE)
    #gIE
    du[8] = alphaIE * (-gIE + pIE)
    #psiEI
    du[9] = alphaEI * (-pEI +  rI)
    #gEI
    du[10] = alphaEI * (-gEI + pEI)
    #psiII
    du[11] = alphaII * (-pII +  rI)
    #gII
    du[12] = alphaII * (-gII + pII)
end
# unused functions  - - - -
#function for model with shunts + gaps
function rhs_function_SHUNTS(du,u,p,t)
    rE=u[1]
    rI=u[2]
    vE=u[3]
    vI=u[4]
    pEE=u[5]
    gEE=u[6]
    pIE=u[7]
    gIE=u[8]
    pEI=u[9]
    gEI=u[10]
    pII=u[11]
    gII=u[12]
    #rE
    du[1] =(1/tauE)*(-gEE*rE -gEI*rE - kappaVE * rE +2 * rE * vE + (deltaE / (tauE*pi)))
    #rI
    du[2] =(1/tauI)*(-gII*rI - gEI*rI- kappaVI * rI + 2 * rI * vI + (deltaI / (tauI*pi)))
    #vE
    du[3] =(1/tauE)*(gEE*(VsynEE - vE) + gEI*(VsynEI - vE) - (tauE^2)*(pi^2) * (rE^2) +  vE^2 + eta_0E)
    #vI
    du[4] =(1/tauI)*(gIE*(VsynIE -vI) + gII*(VsynII - vI) - (pi^2)*(tauI^2) * (rI^2) + vI^2 + eta_0I)
    #psiEE
    du[5] = alphaEE * (-pEE + betaEE * rE)
    #gEE
    du[6] = alphaEE * (-gEE + pEE)
    #psiIE
    du[7] = alphaIE * (-pIE + betaIE * rE)
    #gIE
    du[8] = alphaIE * (-gIE + pIE)
    #psiEI
    du[9] = alphaEI * (-pEI + betaEI * rI)
    #gEI
    du[10] = alphaEI * (-gEI + pEI)
    #psiII
    du[11] = alphaII * (-pII + betaII * rI)
    #gII
    du[12] = alphaII * (-gII + pII)
end


#Functions for Kuromoto type model
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
    du[1] = (1/tauE)*(F(u[1], eta_0E, deltaE, kappaVE) +
            G(u[1], u[4], u[8], VsynEE, VsynEI))
    #zI
    du[2] = (1/tauE)*(F(u[2], eta_0I, deltaI, kappaVI) +
            G(u[2], u[10], u[6], VsynII, VsynIE))
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
