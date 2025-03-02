close all
clear

%% Constants
g = 32.174; % Gravity in ft/s^2
rho = 0.0023769; % Air density at sea level in slugs/ft^3

%% Aircraft Specifications
AR = 2.45; % Aspect ratio
Cl_max = 1.12; % Maximum lift coefficient with high-lift devices
Cd0 = 0.01; % Zero-lift drag coefficient
e = 0.8; % Oswald efficiency factor

%% Derived Parameters
W_S = linspace(0, 100, 500); % Wing loading range in lbs/ft^2
q = 0.5 * rho * (1.467 * 250)^2; % Dynamic pressure at 250 knots (ft/s)

%% Thrust-to-Weight Ratio Calculations
% Takeoff Constraint
Sg = 2000; % Takeoff ground roll in ft
mu = 0.04; % Rolling friction coefficient
T_W_takeoff = (1.21 / (g * rho * Cl_max * Sg)) * W_S + (0.605 / Cl_max) * (Cd0 - mu * Cl_max) + mu;

% Climb Constraint
V_v = 79.016 + 1.2722*(W_S); % Rate of climb in KCAS
V_inf = 399.843; % Climb speed in KCAS
q_climb = 0.5 * rho * V_inf^2;
cd_min = 0.02;
k = 1 / (pi * e * AR);
T_W_climb = (V_v / V_inf) + (q_climb ./ W_S) * cd_min + (k/q_climb)*(W_S);

% Cruise Constraint
V_cruise = 660; % Cruise speed in ft/s (approximately Mach 0.85)
rho = 7.382e-4;% Density of Air at 35k feet
q_cruise = 0.5 * rho * V_cruise^2;
n = 1;
T_W_cruise = q_cruise * (cd_min ./ (W_S) + k*(n/q_cruise)^2 * (W_S));

% Sustained Turn at Mach 0.9 Constraint
n_09M = 4; % Load factor for sustained turn at M=0.8
V_09M = 875.7; % Speed in ft/s (approximately Mach 0.8)
rho = 7.382e-4;% Density of Air at 35k feet
q_09M = 0.5 * rho * V_09M^2;
alpha = 0.31; % Ratio of Air Densities
T_W_turn_09M = (q_09M / alpha) .* ((cd_min ./ W_S) + k * ((n_09M*0.757 / q_09M).^2).*W_S);

% Sustained Turn at Mach 1.2 Constraint
n_12M = 3; % Load factor for sustained turn at M=1.2
V_12M = 1167; % Speed in ft/s (approximately Mach 1.2)
rho = 7.382e-4;% Density of Air at 35k feet
k = 1/(pi * 0.4 * AR);
q_12M = 0.5 * rho * V_12M^2;
T_W_turn_12M = (q_12M / alpha) .* ((cd_min ./ W_S) + k * ((n_12M*0.757 / q_12M).^2).*W_S);

% Mach 1.6 Dash Constraint
n = 2; 
V_16M = 1556; % Speed in ft/s (approximately Mach 1.4)
q_16M = 0.5 * rho * V_16M^2;
T_W_dash_16M = q_16M * (cd_min ./ (W_S) + k*(n/q_16M)^2 * (W_S));

% Landing Constraint
S_LDG = 3500; % Landing distance in feet
h_obst = 50; % Obstacle height in feet
tau = 0.05; % Approach speed parameter (assumed)
C_L_max = 1.2; % Maximum lift coefficient
C_D_LDG = 0.15; % Drag coefficient during landing
mu = 0.4; % Runway friction coefficient
T_gr = 0; % Thrust during ground roll (assumed negligible)
W = 1; % Normalized weight (it cancels out in equations)

% Landing distance equation

%% Plotting

color = parula(7);

figure(1);
hold on;
plot(W_S, T_W_takeoff, "LineWidth", 2,'Color',color(1,:));
plot(W_S, T_W_climb, "LineWidth", 2,'Color',color(2,:));
plot(W_S, T_W_cruise, "LineWidth", 2,'Color',color(3,:));
plot(W_S, T_W_turn_09M, "LineWidth", 2,'Color',color(4,:));
plot(W_S, T_W_turn_12M, "LineWidth", 2,'Color',color(5,:));
plot(W_S, T_W_dash_16M, "LineWidth", 2,'Color',color(6,:));
xline(83.3, "LineWidth", 2,'Color',color(7,:));

ax = gca;
ax.FontSize = 16;
xlabel('Wing Loading (W/S) [lbs/ft^2]','FontSize',18);
ylabel('Thrust-to-Weight Ratio (T/W)','FontSize',18);
title('First Design Constraint Diagram','FontSize',20);
legend('Takeoff', 'Climb', 'Cruise', 'Sustained Turn (M=0.9)', 'Sustained Turn (M=1.2)', 'Dash (M=1.6)', 'Landing','FontSize',14,'location','n');
grid on;
hold off;
axis([0 100 0 2])
