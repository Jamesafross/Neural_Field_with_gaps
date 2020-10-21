function FR(z,tau)
    return (tau/(pi))*((1 .-abs.(z).^2)/(abs.(1 .+z).^2))
end
