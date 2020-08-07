using NLsolve

function f!(F, x)
    F[1] = -kappaV *x[1] .+ 2 * x[1] .* x[2] .+ Delta / (tau * pi)
    F[2] = 0 .- (pi * tau * x[1]) .^ 2 .+ x[2].^ 2 .+ eta_0
end

function SteadyState()
    SS = nlsolve(f!, [ 1.0; 1.2])
    return SS.zero
end
