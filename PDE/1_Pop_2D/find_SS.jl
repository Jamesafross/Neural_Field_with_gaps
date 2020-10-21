using NLsolve
# functions used for finding the homogenous steady state
# of the variables R and V
function f!(F, x)
    F[1] = -kappaV *x[1] .+ 2 * x[1] .* x[2] .+ Delta / (tau * pi) #R
    F[2] = 0 .- (pi * tau * x[1]) .^ 2 .+ x[2].^ 2 .+ eta_0        #V
end
function SteadyState()
    SS = nlsolve(f!, [ 2.0; 1])
    return SS.zero
end
