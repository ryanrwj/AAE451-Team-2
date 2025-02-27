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

% === Sizing Variables ===
W0_S = 1;           % lb/ft^2, wing loading
T_W0 = 1;           % thrust to weight ratio
AR = 1;             % wing aspect ratio
M_max = 1;          % maximum Mach number
M_cruise = 1;       % cruise Mach number
C = 1;              % 1/s, specific fuel consumption
h_cruise = 1;       % ft, cruise altitude
e0 = 1;             % Oswald efficiency factor
L_Dmax = 1;         % maximum lift-to-drag ratio
phi_dot = 1;        % rad/s, turn rate
n_missiles = 4;     % number of missiles
C_D0 = 1;           % zero-lift drag coefficient

% === Mission Requirements ===

% DCA Mission
h_mission = 35000;  % ft, mission altitude
R_cruiseOut = 300;  % nm, at optimum speed/altitude
t_loiter    = 4.0;  % hours, air patrol at 35,000 ft
R_dash      = 100;  % nm, maximum speed dash
x_combat    = 2;    % number of combat turns
R_cruiseBack= 400;  % nm, return to base
t_reserve   = 0.5;  % hours, reserve loiter

% === Initial Guess or Known Values ===
% You might guess an initial gross weight or base empty weight from F-15

W_missile = 356;                    % lb, AIM-120 weight
W_payload = n_missiles*W_missile;   % lb, total weapons
W0_guess = 27853;                   % lb, initial guess for GTOW

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

% For demonstration, let’s define some typical speeds:
V_cruise = M_cruise*a_cruise;

V_cruise  = 450; % knots (placeholder)
V_loiter  = 300; % knots (placeholder)

% Convert knots to nm/hr if needed
% 1 knot = 1 nm/hr, so V_cruise in nm/hr = 450 nm/hr (conveniently the same).

%% 4) Weight Fractions for Each Mission Segment

% --- Segment 1: Taxi + Take-off ---
wf_takeoffFraction = 0.98; % Raymer 6.8

% --- Segment 2: Climb ---
wf_climbFraction = 1.0065 - 0.325*M_cruise; % 6.9

% --- Segment 3: Cruise Out (300 nm) ---
q_cruise = 0.5*rho_cruise*V_cruise^2;   % dynamic pressure at cruise CHECK UNITS
L_Dcruise = ((q_cruise*C_D0/W0_S) + W0_S/(q_cruise*pi*e0*AR))^-1; % 6.13

wf_cruiseOutFraction = exp((-R_cruiseOut*C)/(V_cruise*L_Dcruise)); % 6.11

% --- Segment 3: Loiter (4 hours at 35,000 ft) ---
%   Breguet endurance: Wi/Wf = exp( (E * SFC) / (L/D) )
%   => wf_loiter = Wf/Wi
%
Wi_loiter = 1;
wf_loiter = exp( -(t_loiter * C_T) / (LoverD_loiter));

% --- Segment 4: Dash to target (100 nm, max speed) ---
%   Possibly use a fraction from historical data or Breguet with new speed & L/D
%   We have wf_dashFraction = 0.95 from historical data. 

% --- Segment 5: Combat (2 x 360 turns, missile shots, etc.) ---
%   Could be a small fraction, e.g., 0.97 leftover (wf_combatFraction).

% --- Segment 6: Cruise Back (400 nm) ---
%   Use Breguet again
Wi_cruiseBack = 1;
wf_cruiseBack = exp( -(R_cruiseBack * C_T) / (V_cruise * LoverD_cruise) );

% --- Segment 7: Reserve / Loiter (0.5 hr) ---
Wi_reserve = 1;
wf_reserve = exp( -(t_reserve * C_T) / (LoverD_loiter) );

%% 5) Combine Segment-by-Segment Weight Fractions
% Start with W0, multiply all fractions for each segment to get final weight.

% Let’s define each segment’s fraction in order:
% (Take-off & Climb) -> (Cruise Out) -> (Loiter) -> (Dash) -> (Combat)
% -> (Cruise Back) -> (Reserve)

wf_takeoff_climb = wf_climbFraction;

% Multiply them in sequence:
% Wfinal / W0 = wf_total = wf_takeoff_climb * wf_cruiseOut * wf_loiter 
%                            * wf_dashFraction * wf_combatFraction 
%                            * wf_cruiseBack * wf_reserve

wf_total = wf_takeoff_climb ...
         * wf_cruiseOut ...
         * wf_loiter ...
         * wf_dashFraction ...
         * wf_combatFraction ...
         * wf_cruiseBack ...
         * wf_reserve;

% Then, the total fuel fraction needed is:
% Wf / W0 = 1 - wf_total
fuel_fraction_required = 1 - wf_total;

%% 6) Estimate Gross Weight Based on an Empty Weight Fraction or Iteration
% A common approach:
%   W0 = W_empty + W_crew + W_payload + W_fuel
%   Let We = (some fraction)*W0 from historical data, or use a regression eqn.
%   Then W_fuel = fuel_fraction_required * W0
%
% So W0 = We + W_crew + W_payload + (fuel_fraction_required)*W0
% => W0 - (fuel_fraction_required)*W0 = We + W_crew + W_payload
% => W0 * (1 - fuel_fraction_required) = We + W_crew + W_payload
% => W0 = (We + W_crew + W_payload) / (1 - fuel_fraction_required)
%
% If you know We ~ k*W0, you have an implicit equation that can be solved
% by iteration or using an initial guess. For a first pass, we might 
% approximate We from the F-15 as ~ 0.55 * W0. 
% Then solve iteratively.

% For demonstration, let's do a single iteration approach.

% Step 1: guess We fraction from F-15 data (~0.55)
We_fraction_guess = 0.55;
% Step 2: solve for W0 from the formula above
%   W0 = (We + W_crew + W_payload) / (1 - fuel_fraction_required)
% but We = We_fraction_guess * W0, so
%   W0 = (We_fraction_guess*W0 + W_crew + W_payload) / (1 - fuel_fraction_required)

% Rearrange:
%   W0 - We_fraction_guess*W0/(1 - fuel_fraction_required) = (W_crew + W_payload)/(1 - fuel_fraction_required)
%   W0 * [1 - We_fraction_guess/(1 - fuel_fraction_required)] = (W_crew + W_payload)/(1 - fuel_fraction_required)
%   W0 = [ (W_crew + W_payload)/(1 - fuel_fraction_required ) ] ...
%         / [ 1 - We_fraction_guess/(1 - fuel_fraction_required) ]

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

