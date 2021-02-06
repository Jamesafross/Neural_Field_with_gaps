function rhs_function(du,u,p,t)
    du[1] = (1/tau)*( -kappaV*u[1] + 2*u[1]*u[2] + delta/(tau*pi))
    du[2] = (1/tau)*(kappaS*u[4] - (pi^2)*(u[1]^2)*(tau^2) + u[2]^2  + eta_0)
    du[3] = Alpha*(-u[3] + beta*u[1])
    du[4] = Alpha*(-u[4] + u[3])
end
