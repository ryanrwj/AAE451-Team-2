%% Constants
g = 32.174; % Gravity in ft/s^2
rho = 0.0023769; % Air density at sea level in slugs/ft^3

%% Aircraft Specifications
W_empty = 13865; % Empty weight in lbs
W_max = 29027; % Maximum takeoff weight in lbs
T_max = 15800; % Maximum thrust with afterburner in lbs
S = 196.1; % Wing area in ft^2
AR = 2.45; % Aspect ratio
Cl_max = 1.12; % Maximum lift coefficient with high-lift devices
Cd0 = 0.0172; % Zero-lift drag coefficient
e = 0.8; % Oswald efficiency factor

%% Derived Parameters
W_S = linspace(30, 100, 500); % Wing loading range in lbs/ft^2
q = 0.5 * rho * (1.467 * 250)^2; % Dynamic pressure at 250 knots (ft/s)

%% Thrust-to-Weight Ratio Calculations
% Takeoff Constraint
Sg = 8000; % Takeoff ground roll in ft
mu = 0.04; % Rolling friction coefficient
T_W_takeoff = (1.21 / (g * rho * Cl_max * Sg)) * W_S + (0.605 / Cl_max) * (Cd0 - mu * Cl_max) + mu;

% Climb Constraint
V_v = 79.016 + 1.2722*(W_S); % Rate of climb in KCAS
V_inf = 399.843; % Climb speed in KCAS
q_climb = 0.5 * rho * V_inf^2;
k = 1 / (pi * e * AR);
T_W_climb = (V_v / V_inf) + (q_climb ./ W_S) .* (Cd0 + k * (W_S ./ q_climb).^2);

% Cruise Constraint
V_cruise = 660; % Cruise speed in ft/s (approximately Mach 0.85)
q_cruise = 0.5 * rho * V_cruise^2;
T_W_cruise = (q_cruise ./ W_S) .* Cd0;

% Sustained Turn at Mach 0.8 Constraint
n_08M = 5; % Load factor for sustained turn at M=0.8
V_08M = 880; % Speed in ft/s (approximately Mach 0.8)
q_08M = 0.5 * rho * V_08M^2;
T_W_turn_08M = (q_08M ./ W_S) .* (Cd0 + k * (n_08M^2 * (W_S ./ q_08M).^2));

% Sustained Turn at Mach 1.2 Constraint
n_12M = 4; % Load factor for sustained turn at M=1.2
V_12M = 1320; % Speed in ft/s (approximately Mach 1.2)
q_12M = 0.5 * rho * V_12M^2;
T_W_turn_12M = (q_12M ./ W_S) .* (Cd0 + k * (n_12M^2 * (W_S ./ q_12M).^2));

% Mach 1.4 Dash Constraint
V_14M = 1540; % Speed in ft/s (approximately Mach 1.4)
q_14M = 0.5 * rho * V_14M^2;
T_W_dash_14M = (q_14M ./ W_S) .* Cd0;

%% Plotting
figure;
hold on;
plot(W_S, T_W_takeoff);
plot(W_S, T_W_climb);
plot(W_S, T_W_cruise);
plot(W_S, T_W_turn_08M);
plot(W_S, T_W_turn_12M);
plot(W_S, T_W_dash_14M);

xlabel('Wing Loading (W/S) [lbs/ft^2]');
ylabel('Thrust-to-Weight Ratio (T/W)');
title('First Design Constraint Diagram');
legend('Takeoff', 'Climb', 'Cruise', 'Sustained Turn (M=0.8)', 'Sustained Turn (M=1.2)', 'Dash (M=1.4)');
grid on;
hold off;
