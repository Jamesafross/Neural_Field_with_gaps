
include("firingrate.jl")


function SingleNodeEM(u0,p,T,dt,saveat)
    N = Int(T/saveat)
    StepSave = Int(saveat/dt)
    sol = complex(zeros(length(u0),N))
    u = u0
    c = 1
    numODEs = length(u0)
    du = complex(zeros(size(u0)))

    for i = 1:Int(T/dt)
        @fastmath u .= u .+ dt*(f(u,du,p))

        if mod(i-1,StepSave) == 0
            sol[:,c] .= u[:]
            c += 1
        end
    end
    return sol
end
