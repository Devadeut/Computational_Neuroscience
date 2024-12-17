clear; clc; close all;

%% Parameters
gca = 4.4; gk = 8; gl = 2;
vca = 120; vk = -84; vl = -60;
C = 20; 
Iext = 0;    % No external current
phi = 0.02;
V1 = -1.2; V2 = 18; V3 = 2; V4 = 30;

% Steady-state functions
M_inf = @(V) 0.5*(1 + tanh((V - V1)/V2));
W_inf = @(V) 0.5*(1 + tanh((V - V3)/V4));

% Find equilibrium
equations = @(x) [
    Iext - gca*M_inf(x(1))*(x(1)-vca) - gk*x(2)*(x(1)-vk) - gl*(x(1)-vl);
    W_inf(x(1)) - x(2)
];

x0 = [-60; 0.01];
options = optimoptions('fsolve','Display','off','TolFun',1e-12,'TolX',1e-12);
[x_eq, ~] = fsolve(equations, x0, options);

V_star = x_eq(1);
w_star = x_eq(2);
fprintf('Equilibrium: V* = %.4f mV, w* = %.4f\n', V_star, w_star);

%% ODE function
function dydt = mle_system(~, y, gca, gk, gl, vca, vk, vl, C, phi, Iext, V1, V2, V3, V4)
    V = y(1);
    w = y(2);
    M_inf = 0.5*(1 + tanh((V - V1)/V2));
    W_inf = 0.5*(1 + tanh((V - V3)/V4));
    tau_w = 1./cosh((V - V3)/(2*V4));
    
    dVdt = (Iext - gca*M_inf*(V - vca) - gk*w*(V - vk) - gl*(V - vl))/C;
    dwdt = phi*(W_inf - w)./tau_w;
    dydt = [dVdt; dwdt];
end

%% Simulation parameters
tspan = [0 300]; % ms

% Scan a range of initial voltages
 
% This range assumes equilibrium is around -60.9 mV; adjust as needed.
% If equilibrium is ~-60.9
% so let's pick a better range near the suspected threshold:
% Let's say from -53 mV to -48 mV absolute, not relative to V_star.
V_range = -62:0.001:-60; 

maxV = zeros(size(V_range));
for i = 1:length(V_range)
    Vinit = V_range(i);
    [T,Y] = ode45(@(t,y) mle_system(t,y,gca,gk,gl,vca,vk,vl,C,phi,Iext,V1,V2,V3,V4), tspan, [Vinit; w_star]);
    maxV(i) = max(Y(:,1));

end

% Plot max amplitude vs initial V
figure('Color','w');
plot(V_range, maxV,'LineWidth',0.5);
xlabel('Initial V (mV)');
ylabel('Max V reached (mV)');
title('Max Amplitude of Response vs Initial Depolarization');
grid on;

% From inspection, suppose we see that around -52 mV we get no spike, and by about -51.7 mV we do.
% Let's choose these as examples:

no_spike_Vinit = -60.82; % Subthreshold (adjust based on the plot)
spike_Vinit = -60.6;    % Suprathreshold (adjust based on the plot)

[T_sub, Y_sub] = ode45(@(t,y) mle_system(t,y,gca,gk,gl,vca,vk,vl,C,phi,Iext,V1,V2,V3,V4), tspan, [no_spike_Vinit; w_star]);
[T_spike, Y_spike] = ode45(@(t,y) mle_system(t,y,gca,gk,gl,vca,vk,vl,C,phi,Iext,V1,V2,V3,V4), tspan, [spike_Vinit; w_star]);

% Phase plane plot
figure('Color','w');
plot(Y_sub(:,1), Y_sub(:,2), 'b','LineWidth',1, 'DisplayName','No AP');
hold on; grid on;
plot(Y_spike(:,1), Y_spike(:,2), 'r','LineWidth',1,'DisplayName','AP');
xlabel('V (mV)');
ylabel('w (dimensionless)');
title('Phase-Plane Trajectories: Subthreshold vs Suprathreshold');
legend('Location','best');
hold off;

fprintf('Depolarization happens inbetween range of -60.85 to -60.40');
