function fss(r,p)

    v = (1/(-2*r))*(-beta*r^2 - kappa .*r + delta/pi)
    f1 = (-beta*r^2 - kappa*r + 2*r*v + delta/pi)
    f2 = (r*(v_syn-v) - (pi^2)*(r^2) + v^2 + eta_0)

    return f1 - f2

end

function SS(p,a,b)

f(x) = fss(x,p)
r = find_zero(f,(a,b))

ss = [r, (1/(2*r))*(beta*r^2 + kappa*r - delta/pi), beta*r, beta*r]

return ss
end
