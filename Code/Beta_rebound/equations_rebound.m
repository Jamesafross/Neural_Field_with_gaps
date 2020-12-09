function y = equations_rebound(t,x,delta,eta0,kappaV,kappaS,alpha,tau,alphaD,Omega,driveon,driveoff)

z = x(1);
K = x(2);               % dummy variable, K = (1+alpha_1^(-1)d/dt)U
U = x(3);

% drive  -  (1+alpha_2^(-1)d/dt)^2 U_d = Omega
K_d = x(4);
U_d = x(5);

% Order parameter equations (z is cromplex)
dzdt = 1/tau*((-1i/2)*(z-1)^2 + (1/2)*((z+1)^2).*(-delta + 1i*(eta0 + U_d + kappaV*V(z) + kappaS*U)) - (1/2)*(z^2-1)*(kappaV));

% Synapse equations
dkdt = alpha*(-K + f(z,tau));
dgdt = alpha*(-U + K);

% Synapse equations for drve
dksdt = alphaD*(-K_d + drive(t,driveon,driveoff,Omega));
dgsdt = alphaD*(-U_d + K_d);


y = [dzdt;dkdt;dgdt;dksdt;dgsdt];

end

function y = f(z,tau)
% Firing rate function

y = 1/(tau*pi)*(1-abs(z).^2)./(1+z+conj(z)+abs(z).^2);

end

function y = V(z)
% Compute mean membrane potential V

y = imag((1-conj(z))./(1+conj(z)));

end

function drive = drive(t,driveon,driveoff,strength)

if t>=driveon && t<=driveoff
%% Square pulse
    drive = strength;
    
%% Linearly decreasing drive    
%     drive = (driveoff-t)*strength./(driveoff-driveon);

%% Exponentially decaying pulse
%     drive = strength*exp(-(t-driveon)*rate);

else
    drive = 0;
end

end

