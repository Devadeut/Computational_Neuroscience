clear; clc; close all;

%% Parameters for Morris–Lecar
gca = 4.4; gk = 8; gl = 2;
vca = 120; vk = -84; vl = -60;
C = 20;
Iext = 90; % A current large enough to induce a spike, adjust if needed
V1 = -1.2; V2 = 18; V3 = 2; V4 = 30;

% Define the gating functions
M_inf = @(V) 0.5*(1 + tanh((V - V1)/V2));
W_inf = @(V) 0.5*(1 + tanh((V - V3)/V4));
tau_w = @(V) 1./cosh((V - V3)/(2*V4));

% Simulation parameters
tspan = [0 200]; % ms
V0 = -60; w0 = 0.01;

% First run with phi = 0.02
phi_1 = 0.02;
[~, Y1] = ode45(@(t,y) mle_equations(t,y,phi_1), tspan, [V0; w0]);
V1 = Y1(:,1); w1 = Y1(:,2);

% Second run with phi = 0.04
phi_2 = 0.04;
[~, Y2] = ode45(@(t,y) mle_equations(t,y,phi_2), tspan, [V0; w0]);
V2 = Y2(:,1); w2 = Y2(:,2);

% Plot V vs t for both phi values
figure('Color','w');
subplot(2,1,1);
plot(tspan(1):0.1:tspan(2),0); % just a baseline
hold on; grid on;
[tt, ~] = ode45(@(t,y) mle_equations(t,y,phi_1), tspan, [V0; w0]);
plot(tt, interp1(tt,V1,tt), 'b', 'DisplayName','\phi=0.02');
[tt, ~] = ode45(@(t,y) mle_equations(t,y,phi_2), tspan, [V0; w0]);
plot(tt, interp1(tt,V2,tt), 'r', 'DisplayName','\phi=0.04');
xlabel('time (ms)');
ylabel('V (mV)');
title('Action Potentials for Different \phi');
legend('Location','best');

% Phase plane plot
subplot(2,1,2);
plot(V1, w1, 'b', 'LineWidth', 2, 'DisplayName','\phi=0.02');
hold on; grid on;
plot(V2, w2, 'r', 'LineWidth', 2, 'DisplayName','\phi=0.04');

xlabel('V (mV)');
ylabel('w (dimensionless)');
title('Phase Plane Trajectories');
legend('Location','best');

% Now, as a check, try phi = 0.01
phi_3 = 0.01;
[~, Y3] = ode45(@(t,y) mle_equations(t,y,phi_3), tspan, [V0; w0]);
V3 = Y3(:,1); w3 = Y3(:,2);

figure('Color','w');
plot(V1, w1, 'b', 'LineWidth', 2, 'DisplayName','\phi=0.02');
hold on; grid on;
plot(V2, w2, 'r', 'LineWidth', 2, 'DisplayName','\phi=0.04');
plot(V3, w3, 'g', 'LineWidth', 2, 'DisplayName','\phi=0.01');
xlabel('V (mV)');
ylabel('w');
title('Phase Plane with Different \phi Values');
legend('Location','best');

function dydt = mle_equations(~,y,phi)
    % Unpack variables
    V = y(1);
    w = y(2);
    
    % Global parameters used
    gca = 4.4; gk = 8; gl = 2;
    vca = 120; vk = -84; vl = -60;
    C = 20; Iext = 90;
    V1 = -1.2; V2 = 18; V3 = 2; V4 = 30;
    
    M_inf = @(V) 0.5*(1 + tanh((V - V1)/V2));
    W_inf = @(V) 0.5*(1 + tanh((V - V3)/V4));
    tau_w = @(V) 1./cosh((V - V3)/(2*V4));

    dVdt = (Iext - gca*M_inf(V)*(V - vca) - gk*w*(V - vk) - gl*(V - vl))/C;
    dwdt = phi*(W_inf(V)-w)/tau_w(V);

    dydt = [dVdt; dwdt];
end
