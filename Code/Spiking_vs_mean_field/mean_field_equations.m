function y = mean_field_equations(t,x,delta,eta0,kappaV,kappaS,alpha,tau)

z = x(1);
K = x(2);               % dummy variable, K = (1+alpha_1^(-1)d/dt)g 
G = x(3);

% Order parameter equations (z is complex)
dzdt = 1/tau*((-1i/2)*(z-1)^2 + (1/2)*((z+1)^2).*(-delta + 1i*(eta0 + kappaV*V(z) + kappaS*G)) - (1/2)*(z^2-1)*(kappaV));

% Synapse equations
dkdt = alpha*(-K + f(z,tau));
dgdt = alpha*(-G + K);

y = [dzdt;dkdt;dgdt];

end

function y = f(z,tau)
% Firing rate function f

y = 1/(tau*pi)*(1-abs(z).^2)./(1+z+conj(z)+abs(z).^2);

end

function y = V(z)
% Function to compute mean voltage

y = imag((1-conj(z))./(1+conj(z)));

end

