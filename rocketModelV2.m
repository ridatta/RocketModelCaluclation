function rocketModelV2()
clc; close all; clear;
% Constants
load physicalConstants-SI.mat mu0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Array Parameters (CHANGE THIS)
R0 = 16e-3 / 2; % Array radius, [m]
V = 100 * 10^3; % [m/s]
alph = 0.5; % fraction of mass remaining at t = t0 (end of current pulse)
N = 16; % number of wires
rho_w = 2710; % wire density, [kg/m^3]
t0 = 480e-9; % pulse duration , s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calulate for a given wire diameter
d = 40e-6; % [m]

% % mass in array
rw = 75e-6 / 2; % wire radius [m]
array_mass = N * pi * rw^2 * 1 * rho_w; 

fn = @(t) I(t).^2; % fn to integrate

% Initial mass required

m0 = mu0 / (4 * pi * R0 * V) / (1-alph) * integral(fn,0,t0); % [kg/m] mass per unit length
d_w = getWireDiameter(m0,N,rho_w);  % Initial Wire radius, [m]

% Output
fprintf('m0/l required [x 10^3 kg/m] = %f\n', m0*1e3); 
fprintf('d_w required for %2.0f wires [um] = %f\n', N, d_w * 10^6); 


% show pulse
showCurrentPulse(t0);
saveas(gcf,'current.png');
% % Density distribution
% showDensityDistribution([t0/6, t0/3, t0/2, 4 * t0/6, 5 * t0/6, t0],R0,V);
% saveas(gcf,'denisty.png');


fprintf('num wires for %2.1f um diameter wires = %2.1f\n',d*1e6, getWireNumrequired(m0, rho_w, d));



% % analytical solution for initial mass for I ~ I0 sin(pi/t0 * t)

% % Current pulse parameters - for analytical solution
% I_pk = 1.4e6; % [A]
% t0 = 2 * 240e-9; % [s], Rise time
% m0_anal = mu0 * I_pk^2 / (4 * pi * R0 * V) / (1-alph) * (t0/2); % analytical solution
% d_w_anal = getWireDiameter(m0_anal,N,rho_w);  % Initial Wire radius, [m]
% fprintf('\nAnalytical solution for I(t) = I sin(pi/t0 * t)\n'); 
% fprintf('m0/l (anal) [kg/m] = %f\n', m0_anal); 
% fprintf('d_w [um] = %f\n', d_w_anal * 10^6); 
end


function out = I(t) % Current function  
% Approx. current profile for a underdamped function
% t in [s]

% % % Puffin 
% I_pk = 1 * 10^6; % [A]
% t0 = 2 * 1.5 * 10^-6; % [s], 2 \times rise time
% out = I_pk * sin(pi * t / t0); % [A]
% 
% % Puffin - Exact
% % V0 = 70e3; % [V], Voltage
% % L = 83e-9; % [H], Inductance
% % w = 5.35 * 10^5; % [s^-1]
% % gamma = 3.01 * 10^5; % [s^-1]
% % out = V0 / (L * w) * exp(-gamma*t) .* sin(w * t); % [A]

% Z
Im = 1.4e6; % [A]
tm = 240e-9; % [s], Rise time
out = Im * (sin(pi * (t) / (2 * tm))).^2;

% % MAGPIE
% Im = 1.4e6; % [A]
% tm = 1 * 250e-9; % [s], Rise time
% out = Im * (sin(pi * t / (2 * tm))).^2;

% % MAGPIE sum of sines
% a1 = 0.7987;
% b1 =0.005901;
% c1 = 0.01495;
% a2 = -0.01834;
% b2 = 0.04547 ;
% c2 = -1.248;
% a3 = 0.2095;
% b3 = 0.02026;
% c3 = -3.02;
% 
% f1 = @(x)  a1*sin(b1*x+c1) + a2*sin(b2*x+c2) + a3*sin(b3*x+c3);
% 
% out = 1.0e6 .* f1(t/1e-9);

end

function N = getWireNumrequired(target_mass, rho_w, d_w)
    N = target_mass ./ (rho_w * pi * (d_w/2).^2);
end

function out = getWireDiameter(m0,N,rho_w)
r_w = sqrt(m0 / (pi * N * rho_w)); % Initial Wire radius, [m]
out = 2 * r_w; 
end

function showCurrentPulse(t0)
t = linspace(0,t0,100);
figure
plot(t*1e9,I(t)*1e-6,'k','LineWidth',2.5);
xlabel('time, ns'); ylabel('Current, MA');
formatPlots(); legend off; axis square
end

function showDensityDistribution(t,R0,Vab)
% show density distribution at time v0
mu0 = 4 * pi * 1e-7;
r = linspace(R0,2.5*R0,100); 

figure
for ii =1:length(t)
rho = mu0 ./ (8 * pi^2 * R0 * r * Vab^2) .* (I(t(ii) - (r-R0)/Vab)).^2;
rho((t(ii) - (r-R0)/Vab) < 0 ) = 0;
plot((r-R0)*1e3,rho,'Color',sqclr('r',ii),'LineWidth',2.5,...
    'DisplayName',['t = ' num2str(t(ii)*1e9) ' ns']); hold on; 
end
xlabel('$r-R_0$, mm','Interpreter','latex'); 
ylabel('Mass Density, $kgm^{-3}$','Interpreter','latex');
formatPlots(); axis square
end

