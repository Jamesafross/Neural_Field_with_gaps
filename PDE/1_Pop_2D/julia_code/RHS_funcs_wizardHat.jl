function rhsFun(du, u, p, t)
    df1 = ddR(u[1:XX], u[XX+1:2*XX], u[7*XX+1:8*XX], tau, Delta, kappaS,kappaV, eta_0)
    f1 = F1(u[1:XX], u[XX+1:2*XX],kappaV, Delta, tau)
    A4 = (3/2)*(Dxxyy*u[4*XX+1:5*XX]) .- ((1/v)*f1 .+ (1/v^2)*df1 .- (3/2) * (Dxxyy * u[1:XX] )  )
    du[1:XX] = f1
    du[XX+1:2*XX] = G1(u[1:XX], u[XX+1:2*XX], u[7*XX+1:8*XX], kappaS, tau, eta_0)
    du[2*XX+1:3*XX] = v * (-u[2*XX+1:3*XX] + u[3*XX+1:4*XX])
    du[3*XX+1:4*XX] = v * (-u[3*XX+1:4*XX] + u[4*XX+1:5*XX] + (3/2)*(Dxxyy * u[2*XX+1:3*XX]))
    du[4*XX+1:5*XX] = v * (-u[4*XX+1:5*XX] + u[5*XX+1:6*XX])
    du[5*XX+1:6*XX] = v * (-u[5*XX+1:6*XX] + A4)
    du[6*XX+1:7*XX] = alfa * (-u[6*XX+1:7*XX]+ u[2*XX+1:3*XX])
    du[7*XX+1:8*XX] = alfa * (-u[7*XX+1:8*XX] + u[6*XX+1:7*XX])
end
