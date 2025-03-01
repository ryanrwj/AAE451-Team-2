%% HDI Iterative Parametric Sizing Script
%  Author: Team 2
%  Date:   2/27/25

clear; clc; close all;

%% 1) Define Known/Assumed Inputs

% === General constants ===
g = 32.174;             % ft/s^2, gravitational acceleration
rho_sl = 0.0023769;     % slug/ft^3, sea-level air density
gamma_air = 1.4;        % specific heat ratio of air
R_air = 1717;           % ft-lb/slug-R, gas constant of air
a_cruise = 661;         % Knots
rho_cruise = 0.0019;    % slug/ft^3, at cruise

% === Sizing Variables ===
<<<<<<< Updated upstream
W0_guess = 37500;                   % lb, initial guess for GTOW
wing_area = 300;
W0_S = W0_guess / wing_area;           % lb/ft^2, wing loading
T_W0 = 1.095;           % thrust to weight ratio
AR = 3.2;             % wing aspect ratio
M_max = 2;          % maximum Mach number
M_cruise = 0.9;       % cruise Mach number
C_cruise = 0.8;              % 1/s, specific fuel consumption cruise
C_loiter = 0.35;
h_cruise = 35000;       % ft, cruise altitude
e0 = 0.8;             % Oswald efficiency factor
L_Dmax = 12;         % maximum lift-to-drag ratio
L_Dcruise = 0.0866*L_Dmax;
L_dloiter = L_Dmax;
phi_dot = 18 * pi / 180;        % rad/s, turn rate
=======
W0_S = 1;           % lb/ft^2, wing loading
T_W0 = 1;           % thrust to weight ratio
AR = 1;             % wing aspect ratio
M_max = 1;          % maximum Mach number
V_max = 800;        % kn, maximum speed
M_cruise = 1;       % cruise Mach number
C = 1;              % 1/s, specific fuel consumption
h_cruise = 1;       % ft, cruise altitude
e0 = 1;             % Oswald efficiency factor
L_Dmax = 1;         % maximum lift-to-drag ratio
phi_dot = 1;        % rad/s, turn rate
>>>>>>> Stashed changes
n_missiles = 4;     % number of missiles
C_D0 = 1;           % zero-lift drag coefficient

% === Mission Requirements ===

% DCA Mission
h_mission = 35000;  % ft, mission altitude
R_cruiseOut = 300;  % nm, at optimum speed/altitude
t_loiter    = 4.0;  % hours, air patrol at 35,000 ft
R_dash      = 100;  % nm, maximum speed dash
n_turns    = 2;    % number of combat turns
R_cruiseBack= 400;  % nm, return to base
t_reserve   = 0.5;  % hours, reserve loiter

% === Initial Guess or Known Values ===
% You might guess an initial gross weight or base empty weight from F-15

W_missile = 356;                    % lb, AIM-120 weight
W_payload = n_missiles*W_missile;   % lb, total weapons
W0_guess = 37500;                   % lb, initial guess for GTOW

% For advanced sizing, you can iterate, but we'll do a single pass here.

%% 2) Convert Unit Systems Where Needed
% For instance, Breguet uses consistent units (SI or Imperial).
% Let’s assume you’re working in Imperial units (lb, nm, lb/hr for SFC).
% Adjust as needed.

% Convert distances from nm to miles if you want, or keep them in nm.
% For Breguet in Imperial, you might keep nm if you do the correct SFC and g.
% This example is simplified, so we will keep nm.

% Convert times from hours to seconds if needed:
t_loiter_sec   = t_loiter * 3600; 
t_reserve_sec  = t_reserve * 3600;

%% 3) Breguet Range & Endurance Equations
% Breguet range (for a turbojet/turbofan) in a simplified form (Imperial):
%   R = (V / SFC) * (L/D) * ln(Wi / Wf)
% or an alternate form:
%   Wi / Wf = exp( (R * SFC) / (V * (L/D)) )
%
% Breguet endurance (loiter):
%   E = (1 / SFC) * (L/D) * ln(Wi / Wf)
% or
%   Wi / Wf = exp( (E * SFC) / (L/D) )
%
% You need a speed for cruise, a speed for loiter, or you incorporate them in SFC.

% Calculation for Velocity:
V_cruise = M_cruise*a_cruise;
V_loiter  = 0.7 * V_cruise; % knots 

% Convert knots to nm/hr if needed
% 1 knot = 1 nm/hr, so V_cruise in nm/hr = 450 nm/hr (conveniently the same).

%% 4) Weight Fractions for Each Mission Segment

% --- Segment 1: Taxi + Take-off ---
wf_takeoffFraction = 0.98; % Raymer 6.8

% --- Segment 2: Climb ---
wf_climbFraction = 1.0065 - 0.0325*M_cruise; % 6.9

% --- Segment 3: Cruise Out (300 nm) ---
q_cruise = 0.5*rho_cruise*V_cruise^2;   % dynamic pressure at cruise CHECK UNITS
L_Dcruise = ((q_cruise*C_D0/W0_S) + W0_S/(q_cruise*pi*e0*AR))^-1; % 6.13

wf_cruiseOutFraction = exp((-R_cruiseOut*C)/(V_cruise*L_Dcruise)); % 6.11

% --- Segment 4: Loiter (4 hours at 35,000 ft) ---

wf_loiterFraction = exp((-t_loiter_sec*C)/(L_Dmax));

% --- Segment 5: Dash to target (100 nm, max speed) ---
%   Possibly use a fraction from historical data or Breguet with new speed & L/D
%   We have wf_dashFraction = 0.95 from historical data. 
t_dash = R_dash/V_max;
wf_dashFraction = 1 - C*T_W0*t_dash; % 6.16

% --- Segment 6: Combat (2 x 360 turns, missile shots, etc.) ---
%   Could be a small fraction, e.g., 0.97 leftover (wf_combatFraction).
t_combat = 2*pi*n_turns/phi_dot; % 6.17
wf_combatFraction = 1 - C*T_W0*t_combat; % 6.16

% --- Segment 7: Cruise Back (400 nm) ---

wf_cruiseBackFraction = exp((-R_cruiseBack*C)/(V_cruise*L_Dcruise)); % 6.11

% --- Segment 7: Reserve / Loiter (0.5 hr) ---

wf_reserveFraction = exp((-t_reserve_sec*C)/(L_Dmax));

%% 5) Combine Segment-by-Segment Weight Fractions
% Start with W0, multiply all fractions for each segment to get final weight.

% Let’s define each segment’s fraction in order:
% (Take-off & Climb) -> (Cruise Out) -> (Loiter) -> (Dash) -> (Combat)
% -> (Cruise Back) -> (Reserve)

% Multiply them in sequence:
% Wfinal / W0 = wf_total = wf_takeoff_climb * wf_cruiseOut * wf_loiter 
%                            * wf_dashFraction * wf_combatFraction 
%                            * wf_cruiseBack * wf_reserve

wf_total = wf_takeoffFraction ...
         * wf_climbFraction ...
         * wf_cruiseOutFraction ...
         * wf_loiterFraction ...
         * wf_dashFraction ...
         * wf_combatFraction ...
         * wf_cruiseBackFraction ...
         * wf_reserveFraction;

% Then, the total fuel fraction needed is:
% Wf / W0 = 1 - wf_total
fuel_fraction_required = 1 - wf_total;

%% 6) Estimate Gross Weight Based on an Empty Weight Fraction or Iteration

W0 = W0_guess;
We_fraction = a + b*W0^c1*AR^c2*T_W0^c3*W0_S^c4*M_max^c5; 

num = (W_crew + W_payload) / (1 - fuel_fraction_required);
den = 1 - We_fraction_guess/(1 - fuel_fraction_required);
W0_est = num / den;

% Alternatively, do a direct iteration. For a quick example, let's do a 
% simpler approach: We'll guess W0 and see if it is consistent. 
% We'll just demonstrate the direct formula approach:

disp('---------------------');
disp('Preliminary Sizing Results:');
disp(['Total Fuel Fraction Required = ', num2str(fuel_fraction_required)]);
disp(['Estimated Gross Takeoff Weight (W0) = ', num2str(W0_est), ' lb']);

% Now you can compute actual empty weight:
We_est = We_fraction_guess * W0_est;
disp(['Estimated Empty Weight (We) = ', num2str(We_est), ' lb']);

% Fuel weight:
Wf_est = fuel_fraction_required * W0_est;
disp(['Estimated Fuel Weight (Wf) = ', num2str(Wf_est), ' lb']);

disp('---------------------');

%% 7) Check for Convergence or Iterate (Optional)
% In practice, you might:
%  1) Use W0_est to re-compute We_fraction (from a regression eqn or advanced method).
%  2) Recompute fuel fraction with updated aerodynamic or propulsion data.
%  3) Iterate until converged.

% For demonstration, we stop here.

