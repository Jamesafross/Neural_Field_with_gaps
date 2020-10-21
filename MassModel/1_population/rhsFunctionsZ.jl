function V(z)
    return imag((1 - conj(z)) / (1 + conj(z)))
end

function FR(z)
    return (1 / pi) * real((1 - conj(z)) / (1 + conj(z)))
end

function F(z, eta0, delta, kappa)
    return -0.5 * im * ((z - 1)^2) +
           0.5 * ((z + 1)^2) * (-delta + im * eta0 + im * kappa * V(z)) -
           0.5 * (((z^2) - 1)) * kappa
end

function G(z, g1, v_syn1)
    return 0.5 * im * ((z + 1)^2) * v_syn1 * g1 - 0.5 * ((z^2) - 1) * g1

end


function rhs_functionZ(du, u, p, t)
    du[1] = F(u[1], eta_0, delta, kappa) + G(u[1], u[2], v_syn)
    du[2] = Alpha * (-u[2] + beta * FR(u[1]))
    du[3] = Alpha * (-u[3] +  u[2])

end
