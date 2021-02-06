function rhsFun(du, u, p, t)
    # f1 = ∂R/∂t
    # df1 = ∂²R/∂t²
    R .= u[1:XX]
    V .= u[XX+1:2*XX]
    psi .= u[2*XX+1:3*XX]
    A1 .= u[3*XX+1:4*XX]
    A2 .= u[4*XX+1:5*XX]
    A3 .= u[5*XX+1:6*XX]
    p .= u[6*XX+1:7*XX]
    g .= u[7*XX+1:8*XX]

    d2Rdt2 = d2Rdt2Func(R, V, g, tau, Delta, kappaS,kappaV, eta_0)
    dRdt = dRdtFunc(R, V,kappaV, Delta, tau)
    A4 = (3/2)*(Dxxyy*A2) .- ((1/v)*f1 .+ (1/v^2)*df1 .- (3/2) * (Dxxyy * R )  )
    du[1:XX] = f1
    du[XX+1:2*XX] = dVdtFunc(R, V, g, kappaS, tau, eta_0)
    du[2*XX+1:3*XX] = v * (-psi + A1)
    du[3*XX+1:4*XX] = v * (-A1 + A2 + (3/2)*(Dxxyy * psi))
    du[4*XX+1:5*XX] = v * (-A2 + A3)
    du[5*XX+1:6*XX] = v * (-A3 + A4)
    du[6*XX+1:7*XX] = alfa * (-p+ psi)
    du[7*XX+1:8*XX] = alfa * (-g + p)
end
