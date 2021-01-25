clear all
close all

% Model parameters
delta = 0.5;
alpha = 0.1;
kappaV = 1;     % 0.5 for low gap-junction coupling, 1 for intermediate and 1.5 for high
kappaS = 1;
eta0 = 1;
tau = 15;
alphaD = 0.2;
Omega = 3;
driveon = 2000;
driveoff = 2400;

% Simulation parameters
tstart = 0;
tfinal = 5000;  
refine = 4;
options = odeset('OutputSel',1,'Refine',refine);

% Setup initial conditions
z0 = -0.1475 + 0.0165i ; 
U0 = [z0,0.0268,0.0230,0,0];

% Solve mean field equations
sol =  ode23(@equations_rebound,[tstart,tfinal],U0,options,delta,eta0,kappaV,kappaS,alpha,tau,alphaD,Omega,driveon,driveoff);

% Evaluate solution on given interval
t = linspace(0,tfinal,10000);
y = deval(sol,t);

% Define order parameter Z and synaptic current kappaS*U
Z = y(1,:);
current =  kappaS*y(3,:);

% Transform to QIF framework and compute firing rate
W = (1-conj(Z))./(1+conj(Z));
rate = real(W);

% Compute time frequency spectrogram of synaptic current
dt = (t(2)-t(1))/1000;
Fs = 1/(dt);                    % sampling rate in Hz
window = 2500;
overlap = 2490;
[~,F,T,P] = spectrogram(current,window,overlap,10:0.01:35,Fs,'yaxis');
power = 1000*2*abs(P/(window*dt));

% Plot time frequency spectrogram of current
figure(1)
surf(T,F,power,'EdgeColor','none')
view(2)
axis([1.4 4 10 35])
colorbar
set(gca,'linewidth',1.5,'fontsize',24,'fontname','Times')
xlabel('Time (s)','interpreter','latex','FontSize', 28)
ylabel('Frequency (Hz)','interpreter','latex','FontSize', 28)
set(gca,'xtick',[1 2 3 4])
set(gca,'xticklabel',{-1000 0 1000 2000 3000})

% Plot synchrony |Z|
figure
set(gcf,'units','centimeters','position',[1,1,20,10]);
plot(t,abs(Z),'b','linewidth',2)
axis([1400 4000 0 1])
set(gca,'linewidth',1.5,'fontsize',24,'fontname','Times')
xlabel('Time (ms)','interpreter','latex','FontSize', 28)
ylabel('$|Z|$','interpreter','latex','FontSize', 28)
set(gca,'xtick',[1000 2000 3000 4000])
set(gca,'xticklabel',{-1000 0 1000 2000 3000})

% Plot synaptic current
figure
set(gcf,'units','centimeters','position',[1,1,20,10]);
plot(t,current,'r','linewidth',2)
axis([1400 4000 0.01 0.055])
set(gca,'linewidth',1.5,'fontsize',24,'fontname','Times')
xlabel('Time (ms)','interpreter','latex','FontSize', 28)
ylabel('Synaptic current','interpreter','latex','FontSize', 28)
set(gca,'xtick',[1000 2000 3000 4000])
set(gca,'xticklabel',{-1000 0 1000 2000 3000})


