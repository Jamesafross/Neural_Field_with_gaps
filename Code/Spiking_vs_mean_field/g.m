function g = g(centre, half_width, dim)
% Randomly draws drives eta_i from Lorentzian distribution

pd = makedist('tLocationScale','mu',centre,'sigma',half_width,'nu',1);

g = random(pd,dim,1);

end

