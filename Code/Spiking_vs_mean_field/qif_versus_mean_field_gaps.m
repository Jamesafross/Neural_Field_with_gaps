clear all
close all

% Model parameters
delta = 0.5;
alpha = 0.5;
kappaV = 0.8;
kappaS = 1;
eta0 = 2;
tau = 16;
vth = 1000;
vreset = -1000;

% Simulation parameters
dim = 100;                  % number of neurons
tstart = 0;
tfinal = 700;
dt = 1/tau;
refine = 4;

% 'spike.m' defines the threshold events, i.e. v crosses vth. MATLAB's
% event detector runs the ode solver until such an event occurs. You then
% need to reset v and start the ode solver again. You
options = odeset('Events',@spike_gaps,'OutputSel',1,'Refine',refine);

% Setup vectors for spike times etc.
neuron_spike = [];
spike_times = [];
t_all = [];
y_all = [];

% Setup initial conditions
theta = [-0.9*pi*rand(dim/2,1);0.9*pi*rand(dim/2,1)];
y0 = tan(theta/2);
y=[y0',0,0];

% Select background drives from Lorentzian with mean eta0 and width delta
eta = g(eta0, delta, dim);
disp(median(eta))

t = tstart;

for i = 1:dim*20
    if t(end)<tfinal-dt
    % Run until there's been, on average, 50 spikes per neuron or the time
    % reaches tfinal

        % Solve until the first terminal event.
        [t,y,te,ye,ie] = ode23(@qif_gaps_equations,tstart:dt:tfinal,y(end,:),options,eta,dim,kappaV,kappaS,alpha,tau,vth);
        
        % Keep track of neurons which spiked and spike time
        neuron_spike = [neuron_spike; ie];    % Neuron number
        spike_times = [spike_times; te];    % Spike time
        
        % Synpatic input to all neurons. Neurons all-to-all connected so
        % no need to track synapses separately, they're all the same.
        y(end,dim+1)=y(end,dim+1)+alpha/dim;
        
        % Use vi=0 for neuron that spiked, average between vth and vreset
        yend = y(end,:);
        yend(ie) = 0;
        
        % Keep track of time and variables
        t_all = [t_all;t(2:end-1);t(end)];
        y_all = [y_all;y(2:end-1,:);yend];
        
        % Reset voltage of neuron that spiked
        y(end,ie) = vreset;
        
        % Time at which to start at on next loop
        tstart = t(end);

    end
end

% Transform to theta neuron framework
theta = 2*atan(y_all(:,1:dim));
z_spike = sum(exp(1i*theta),2)/dim;
z0 = sum(exp(1i*theta(1,:)))/dim;

% Solve mean field equations
options = odeset('OutputSel',1,'Refine',refine);
[t,z] =  ode23(@mean_field_equations,0:dt:tfinal,[z0;0;0],options,delta,eta0,kappaV,kappaS,alpha,tau);

% Define order parameter for mean field solution
z_mf = z(:,1);

% Compute firing rate R
R_spike = 1/(tau*pi)*(1-abs(z_spike).^2)./(1+z_spike+conj(z_spike)+abs(z_spike).^2);
R_mf = 1/(tau*pi)*(1-abs(z_mf.^2)./(1+z_mf+conj(z_mf)+abs(z_mf).^2));

% Create raster plot
figure(1)
plot(spike_times,neuron_spike,'.')

% Plot mean membrane potential V
figure(2)
hold on
plot(t_all,mean(y_all(:,1:dim),2))
plot(t,imag((1-conj(z_mf))./(1+conj(z_mf))))

% Plot firing rate R
figure
hold on
plot(t_all,R_spike)
plot(t,R_mf)

% Plot synchrony |Z|
figure
hold on
plot(t_all,abs(z_spike))
plot(t,abs(z_mf))



