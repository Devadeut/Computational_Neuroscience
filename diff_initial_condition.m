clear; clc; close all;

%% Parameters
gca = 4.4; gk = 8.0; gl = 2.0;
vca = 120; vk = -84; vl = -60;
phi = 0.02;
V1 = -1.2; V2 = 18; V3 = 2; V4 = 30;
C = 20;

M_inf = @(V) 0.5*(1 + tanh((V - V1)/V2));
W_inf = @(V) 0.5*(1 + tanh((V - V3)/V4));

% ODE function for MLE
function dydt = mle_equations(~,y,gca,gk,gl,vca,vk,vl,C,phi,Iext,V1,V2,V3,V4)
    V = y(1);
    w = y(2);
    M = 0.5*(1 + tanh((V - V1)/V2));
    W = 0.5*(1 + tanh((V - V3)/V4));
    tau_w = 1./cosh((V - V3)/(2*V4));

    dVdt = (Iext - gca*M*(V - vca) - gk*w*(V - vk) - gl*(V - vl))/C;
    dwdt = phi*(W - w)./tau_w;
    dydt = [dVdt; dwdt];
end

%% 1) Equilibrium for Iext=0
Iext_0 = 0;
eq_equations_0 = @(x) [
    Iext_0 - gca*M_inf(x(1))*(x(1)-vca) - gk*x(2)*(x(1)-vk) - gl*(x(1)-vl);
    W_inf(x(1)) - x(2)
];

x0_guess = [-60; 0.01];
options = optimoptions('fsolve','Display','off','TolFun',1e-12,'TolX',1e-12);
[x_eq_0, ~] = fsolve(eq_equations_0, x0_guess, options);
V_star_0 = x_eq_0(1);
w_star_0 = x_eq_0(2);

%% 2) Equilibrium for Iext=86
Iext_86 = 86;
eq_equations_86 = @(x) [
    Iext_86 - gca*M_inf(x(1))*(x(1)-vca) - gk*x(2)*(x(1)-vk) - gl*(x(1)-vl);
    W_inf(x(1)) - x(2)
];

[x_eq_86, ~] = fsolve(eq_equations_86, x0_guess, options);
V_star_86 = x_eq_86(1);
w_star_86 = x_eq_86(2);

fprintf('Equilibrium at Iext=0: V* = %.4f mV, w* = %.4f\n', V_star_0, w_star_0);
fprintf('Equilibrium at Iext=86: V* = %.4f mV, w* = %.4f\n', V_star_86, w_star_86);

%% Simulation setup
tspan = [0 300]; % 300 ms to see full response

% 1) Start from equilibrium at Iext=0, now run with Iext=86
y0_1 = [V_star_0; w_star_0]; 
[~, Y1] = ode45(@(t,y) mle_equations(t,y,gca,gk,gl,vca,vk,vl,C,phi,Iext_86,V1,V2,V3,V4), tspan, y0_1);

% 2) Start from equilibrium at Iext=86, run with Iext=86
y0_2 = [V_star_86; w_star_86];
[~, Y2] = ode45(@(t,y) mle_equations(t,y,gca,gk,gl,vca,vk,vl,C,phi,Iext_86,V1,V2,V3,V4), tspan, y0_2);

% 3) Start from off-equilibrium at (-27.9, 0.17) with Iext=86
y0_3 = [-27.9; 0.17];
[~, Y3] = ode45(@(t,y) mle_equations(t,y,gca,gk,gl,vca,vk,vl,C,phi,Iext_86,V1,V2,V3,V4), tspan, y0_3);

%% Phase-plane plot
figure('Color','w');
plot(Y1(:,1), Y1(:,2), 'b','LineWidth',1, 'DisplayName','From eq at I_{ext}=0');
hold on; grid on;
plot(Y3(:,1), Y3(:,2), 'g','LineWidth',1, 'DisplayName','From (-27.9,0.17), I_{ext}=86');
plot(Y2(:,1), Y2(:,2), 'r','LineWidth',1, 'DisplayName','From eq at I_{ext}=86');
xlabel('V (mV)');
ylabel('w (dimensionless)');
title('Phase Plane: Three Initial Conditions, I_{ext}=86 μA/cm^2');
legend('Location','best');

figure;
plot(Y2(:,1), Y2(:,2), 'r','LineWidth',1, 'DisplayName','From eq at I_{ext}=86');
xlabel('V (mV)');
ylabel('w (dimensionless)');
title('Phase Plane:Initial Conditions, I_{ext}=86 μA/cm^2');
legend('Location','best');
