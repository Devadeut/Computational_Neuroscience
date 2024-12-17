%% Clear workspace
clear; clc; close all;

%% Parameters
gca = 4.4; gk = 8.0; gl = 2.0;
vca = 120; vk = -84; vl = -60;
phi = 0.02;
V1 = -1.2; V2 = 18;
V3 = 2;    V4 = 30;
C = 20;
Iext = 0;

%% Steady-state functions
M_inf = @(V) 0.5*(1 + tanh((V - V1)/V2));
W_inf = @(V) 0.5*(1 + tanh((V - V3)/V4));

%% Equations for equilibrium
equations = @(x) [
    Iext - gca*M_inf(x(1))*(x(1)-vca) - gk*x(2)*(x(1)-vk) - gl*(x(1)-vl);
    W_inf(x(1)) - x(2)
];

%% Solve equilibrium
x0 = [-60; 0.01];
options = optimoptions('fsolve','Display','off','TolFun',1e-12,'TolX',1e-12);
[x_equilibrium, ~] = fsolve(equations, x0, options);
V_star = x_equilibrium(1);
w_star = x_equilibrium(2);

fprintf('Equilibrium: V* = %.4f mV, w* = %.4f\n', V_star, w_star);

%% Compute nullclines
V_values = linspace(-80, 50, 200);
w_nullcline = W_inf(V_values); % w-nullcline

% V-nullcline: solve 0 = -gca*M_inf(V)*(V - vca) - gk*w*(V - vk) - gl*(V - vl)
% w = [ -gca*M_inf(V)*(V - vca) - gl*(V - vl) ] / [gk*(V - vk)]
w_nullcline_V = zeros(size(V_values));
for i = 1:length(V_values)
    Vtemp = V_values(i);
    numerator = -gca*M_inf(Vtemp)*(Vtemp - vca) - gl*(Vtemp - vl);
    denominator = gk*(Vtemp - vk);
    w_nullcline_V(i) = numerator/denominator;
end

%% Compute vector field for quiver
V_grid = linspace(-80, 50, 20);
w_grid = linspace(0, 1, 20);
[Vg, Wg] = meshgrid(V_grid, w_grid);

dVdt = zeros(size(Vg));
dwdt = zeros(size(Wg));

for i = 1:numel(Vg)
    % dV/dt
    dVdt(i) = (Iext - gca*M_inf(Vg(i))*(Vg(i)-vca) - gk*Wg(i)*(Vg(i)-vk) - gl*(Vg(i)-vl))/C;
    % dw/dt
    tau_w = 1./cosh((Vg(i)-V3)/(2*V4));
    dwdt(i) = phi*(W_inf(Vg(i)) - Wg(i))./tau_w;
end

%% Scaling the w-axis
scaleFactor = 150;
% Scale all w-related variables for plotting
W_plot = Wg * scaleFactor;
w_nullcline_plot = w_nullcline * scaleFactor;
w_nullcline_V_plot = w_nullcline_V * scaleFactor;
w_star_plot = w_star * scaleFactor;
dwdt_plot = dwdt * scaleFactor; % scale the derivative as well

%% Plot
figure('Color','w');

% Plot nullclines
plot(V_values, w_nullcline_plot, 'r', 'LineWidth', 2, 'DisplayName','w-nullcline (w=W_\infty(V))');
hold on;
plot(V_values, w_nullcline_V_plot, 'b', 'LineWidth', 2, 'DisplayName','V-nullcline');

% Plot equilibrium point
plot(V_star, w_star_plot, 'ko', 'MarkerSize', 10, 'MarkerFaceColor','k', 'DisplayName','Equilibrium');

% Annotate equilibrium values on the plot
% We'll place the text slightly to the right and above the equilibrium point
text(V_star+2, w_star_plot+0.05*scaleFactor, ...
    sprintf('V* = %.2f mV\nw* = %.4f', V_star, w_star), ...
    'FontSize', 8, 'BackgroundColor','white', 'EdgeColor','black');

% Quiver plot (vector field)
% Use the scaled w-grid for plotting:
quiver(Vg, W_plot, dVdt, dwdt_plot, 'k');

xlabel('V (mV)');
ylabel('w (x100, dimensionless)');
title('Phase Plane Plot of Morris–Lecar Model (w scaled by 150)');
legend('Location','best');
grid on; hold off;
