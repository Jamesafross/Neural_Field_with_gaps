

function rhsFun(du, u, p, t)
    R = u[1:X]
    V = u[X+1:2*X]
    psi = u[2*X+1:3*X]
    A1 = u[3*X+1:4*X]
    A2 = u[4*X+1:5*X]
    A3 = u[5*X+1:6*X]
    p = u[6*X+1:7*X]
    g = u[7*X+1:8*X]
    df1 = ddR(R, V, g, tau, Delta, kappaV, eta_0)
    ddf1 = dddR(R, V, g, kappaV, tau, Delta, eta_0)
    f1 = F1(R, V, kappaV, Delta, tau)


    A4 = Dxx*A2 - 2 * ((1 / v) * f1 .+
                 (2 / v^2) * df1 .+
                 (1 / v^3) * ddf1 .- 2 * Dxx * R .-(1 / v) * Dxx * f1)
    du[1:X]= f1
    du[X+1:2*X] = G1(R, V, g, tau, eta_0)
    du[2*X+1:3*X] = v * (-psi + A1)
    du[3*X+1:4*X] = v * (-A1 + A2 + Dxx * psi)
    du[4*X+1:5*X] = v * (-A2 + A3)
    du[5*X+1:6*X] = v * (-A3 + A4)
    du[6*X+1:7*X] = alfa * (-p + kappaS*psi)
    du[7*X+1:8*X] = alfa * (-g + p)
end
