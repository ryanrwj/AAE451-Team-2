close all;
clear;

%% Constants
g = 32.174; % Gravity in ft/s^2
rho = 0.0023769; % Air density at sea level in slugs/ft^3
rho_icy = 2.938e-3;

%% Aircraft Specifications
AR = 2.9; % Aspect ratio
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
mu_icy = 0.02;
T_W_takeoff = (1.21 / (g * rho * Cl_max * Sg)) * W_S + (0.605 / Cl_max) * (Cd0 - mu * Cl_max) + mu;
T_W_takeoff_icy = (1.21 / (g * rho_icy * Cl_max * Sg)) * W_S + (0.605 / Cl_max) * (Cd0 - mu_icy * Cl_max) + mu_icy;

% Climb Constraint
rho = 1.496e-3;
V_v = 79.016 + 1.2722 * (W_S); % Rate of climb in KCAS
V_inf = 449.882617; % Climb speed in KCAS
q_climb = 0.5 * rho * V_inf^2;
cd_min = 0.02;
k = 1 / (pi * e * AR);
T_W_climb = (V_v / V_inf) + (q_climb ./ W_S) * cd_min + (k / q_climb) * (W_S);

% Cruise Constraint
V_cruise = 660; % Cruise speed in ft/s (Mach 0.85)
rho = 7.382e-4; % Density at 35k feet
q_cruise = 0.5 * rho * V_cruise^2;
T_W_cruise = q_cruise * (cd_min ./ W_S + k * (1 / q_cruise)^2 .* W_S);

% Sustained Turn Constraints
rho = 7.382e-4; % Air density at 35k feet
alpha = 0.31; % Ratio of Air Densities
n_values = [3, 2]; % Load factors for Mach 0.9 and 1.2
M_values = [0.9, 1.2];
V_values = [875.7, 1167]; % Speed in ft/s
T_W_turn = zeros(length(n_values), length(W_S));

for i = 1:length(n_values)
    q_turn = 0.5 * rho * V_values(i)^2;
    T_W_turn(i, :) = (q_turn / alpha) .* ((cd_min ./ W_S) + k * ((n_values(i) * 0.757 / q_turn).^2) .* W_S);
end

% Mach 1.6 Dash Constraint
V_16M = 1556; % Speed in ft/s (Mach 1.6)
q_16M = 0.5 * rho * V_16M^2;
T_W_dash_16M = q_16M * (cd_min ./ W_S + k * (2 / q_16M)^2 .* W_S);

% SEP Constraints (1-g and 5-g at Military & Max Thrust)
altitudes = [0, 15000]; % ft
rhos = [0.002377, 0.00165]; % Air densities
a_speeds = [1116, 1062]; % Speed of sound (ft/s)
Ps_values = [200, 50; 700, 400; 300, 50]; % SEP for 1-g & 5-g

T_W_SEP = zeros(size(Ps_values, 1), length(W_S), length(altitudes));

for j = 1:length(altitudes)
    rho = rhos(j);
    a = a_speeds(j);
    V = 0.9 * a; % Velocity at Mach 0.9
    q_SEP = 0.5 * rho * V^2;
    
    for i = 1:size(Ps_values, 1)
        Ps = Ps_values(i, j);
        T_W_SEP(i, :, j) = (Ps / V) + (q_SEP ./ W_S) .* (Cd0 + (W_S.^2) ./ (pi * e * AR * q_SEP^2));
    end
end

%% Plotting
figure(1);
hold on;
color = parula(12);

plot(W_S, T_W_takeoff, 'LineWidth', 2, 'Color', color(1, :));
plot(W_S, T_W_climb, 'LineWidth', 2, 'Color', color(3, :));
plot(W_S, T_W_cruise, 'LineWidth', 2, 'Color', color(4, :));
plot(W_S, T_W_turn(1, :), 'LineWidth', 2, 'Color', color(5, :));
plot(W_S, T_W_turn(2, :), 'LineWidth', 2, 'Color', color(6, :));
plot(W_S, T_W_dash_16M, 'LineWidth', 2, 'Color', color(7, :));

% Plot SEP constraints
for j = 1:length(altitudes)
    for i = 1:size(Ps_values, 1)
        plot(W_S, T_W_SEP(i, :, j), '--', 'LineWidth', 2, 'Color', color(8 + i, :));
    end
end

xline(83.3, 'LineWidth', 2, 'Color', color(10, :));
historic = [0.9811816192560177, 77.76088564261696; 0.6354485776805252, 75.90603770995509; 1.1999999999999997, 70.17375856784514; 1.0424507658643327, 71.47203013136726];
scatter(historic(:, 2), historic(:, 1));
scatter(77, 0.9, 'filled', 'MarkerFaceColor', 'k');

ax = gca;
ax.FontSize = 16;
xlabel('Wing Loading (W/S) [lbs/ft^2]', 'FontSize', 18);
ylabel('Thrust-to-Weight Ratio (T/W)', 'FontSize', 18);
title('First Design Constraint Diagram', 'FontSize', 20);

legend({'Takeoff', 'Climb', 'Cruise', ...
        'Sustained Turn (M=0.9)', 'Sustained Turn (M=1.2)', 'Dash (M=1.6)', ...
        '1-g SEP (SL)', '1-g SEP (15k ft)', '5-g SEP (SL)', '5-g SEP (15k ft)', ...
        'Landing'}, ...
        'FontSize', 10, 'Location', 'eastoutside');

grid on;
hold off;
axis([0 100 0 2]);
