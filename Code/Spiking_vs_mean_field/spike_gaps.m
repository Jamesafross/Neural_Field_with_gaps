function [value,isterminal,direction] = spike_gaps(t,y,eta,dim,kappaV,kappaS,alpha,tau,vth)
% Need to include all parameter that are input to ode23, even though
% they're not used here

V = y(1:dim);

value = V-vth ;             % Detects when value = 0
isterminal = ones(dim,1);   % Stop the integration
direction = ones(dim,1);


end

