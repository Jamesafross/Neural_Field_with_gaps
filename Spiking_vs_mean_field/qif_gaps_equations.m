function dydt = qif_gaps_equations(t,y,eta,dim,kappaV,kappaS,alpha,tau,vth)

V = y(1:dim);
h = y(dim+1);
g = y(dim+2);

V_sum = sum(V)/dim;

% Voltage equation
dV = (V.^2 + eta + kappaV*(V_sum-V) + kappaS.*g)/tau;

% Synaptic equations
% Neurons all-to-all connected so no need to track synapses separately, 
% as they're all the same
dh = alpha*(-h);
dg = alpha*(h-g);

    
dydt = [ dV; dh; dg];

end
