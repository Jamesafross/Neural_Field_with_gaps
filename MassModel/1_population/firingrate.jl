function FR(z,tau)
    return (1/(pi*tau))*((1 .-abs.(z).^2)/(abs.(1 .+z).^2))
end
