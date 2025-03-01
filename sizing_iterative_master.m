%% HDI Iterative Parametric Sizing Script
%  Author: Team 2
%  Date:   2/27/25

clear; clc; close all;

%% 1) Define Known/Assumed Inputs

% === General constants ===
g = 32.174;             % ft/s^2, gravitational acceleration
rho_sl = 0.0023769;     % slug/ft^3, sea-level air density
rho_cruise = 0.0019;    % slug/ft^3, cruise air density
a_cruise = 661;         % kn, cruise speed of sound
gamma_air = 1.4;        % specific heat ratio of air
R_air = 1717;           % ft-lb/slug-R, gas constant of air

% === Sizing Variables ===
W0_S = 92;           % lb/ft^2, wing loading
W0_guess = 37500;                   % lb, initial guess for GTOW
wing_area = 300;
T_W0 = 1.095;           % thrust to weight ratio
AR = 3.2;             % wing aspect ratio
M_max = 2;          % maximum Mach number
M_cruise = 0.9;       % cruise Mach number
C_cruise = 0.8;              % 1/h, cruise specific fuel consumption
C_loiter = 0.7;
C_dash = 1.4;
h_cruise = 35000;       % ft, cruise altitude
e0 = 0.8;             % Oswald efficiency factor
L_Dmax = 12;         % maximum lift-to-drag ratio
phi_dot = 18 * pi / 180;        % rad/s, turn rate
n_missiles = 4;             % number of missiles
C_D0 = 0.0016;          % zero-lift drag coefficient

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

% For advanced sizing, you can iterate, but we'll do a single pass here.

%% 2) Convert Unit Systems Where Needed
% For instance, Breguet uses consistent units (SI or Imperial).
% Let’s assume you’re working in Imperial units (lb, nm, lb/hr for SFC).
% Adjust as needed.

% Convert distances from nm to miles if you want, or keep them in nm.
% For Breguet in Imperial, you might keep nm if you do the correct SFC and g.
% This example is simplified, so we will keep nm.


%% 3) Velocity calculations

% Calculation for Velocity:
V_cruise = M_cruise*a_cruise;
V_max = M_max*a_cruise;
V_loiter  = 0.7 * V_cruise; % knots 

%% 4) Weight Fractions for Each Mission Segment

% --- Segment 1: Taxi + Take-off ---
wf_takeoffFraction = 0.98; % Raymer 6.8

% --- Segment 2: Climb ---
wf_climbFraction = 1.0065 - 0.0325*M_cruise; % 6.9

% --- Segment 3: Cruise Out (300 nm) ---
q_cruise = 0.5*rho_cruise*V_cruise^2;   % dynamic pressure at cruise CHECK UNITS
L_Dcruise = ((q_cruise*C_D0/W0_S) + W0_S/(q_cruise*pi*e0*AR))^-1; % 6.13

wf_cruiseOutFraction = exp((-R_cruiseOut*C_cruise)/(V_cruise*L_Dcruise)); % 6.11

% --- Segment 4: Loiter (4 hours at 35,000 ft) ---

wf_loiterFraction = exp((-t_loiter*C_loiter)/(L_Dmax));

% --- Segment 5: Dash to target (100 nm, max speed) ---
%   Possibly use a fraction from historical data or Breguet with new speed & L/D
%   We have wf_dashFraction = 0.95 from historical data. 
t_dash = R_dash/V_max;
wf_dashFraction = 1 - C_dash*T_W0*t_dash; % 6.16

% --- Segment 6: Combat (2 x 360 turns, missile shots, etc.) ---
%   Could be a small fraction, e.g., 0.97 leftover (wf_combatFraction).
t_combat = 2*pi*n_turns/phi_dot/3600; % 6.17
wf_combatFraction = 1 - C_dash*T_W0*t_combat; % 6.16

% --- Segment 7: Cruise Back (400 nm) ---

wf_cruiseBackFraction = exp((-R_cruiseBack*C_cruise)/(V_cruise*L_Dcruise)); % 6.11

% --- Segment 7: Reserve / Loiter (0.5 hr) ---

wf_reserveFraction = exp((-t_reserve*C_loiter)/(L_Dmax));

%% 5) Estimate Gross Weight Based on an Empty Weight Fraction or Iteration

W0 = W0_guess;

% polynomial parameters
a = -0.02;
b = 2.16;
c1 = -0.10;
c2 = 0.20;
c3 = 0.04;
c4 = -0.10;
c5 = 0.08;

    We_fraction = a + b*W0^c1*AR^c2*T_W0^c3*W0_S^c4;
    Wf_total = wf_takeoffFraction ...
         * wf_climbFraction ...
         * wf_cruiseOutFraction ...
         * wf_loiterFraction ...
         * wf_dashFraction ...
         * wf_combatFraction ...
         * wf_cruiseBackFraction ...
         * wf_reserveFraction;
    W0 = W_payload/(1 - Wf_total - We_fraction);

%% 6) Display results

fuel_fraction = 1 - Wf_total;
disp('---------------------');
disp('Preliminary Sizing Results:');
disp(['Total Fuel Fraction Required = ', num2str(fuel_fraction)]);
disp(['Estimated Gross Takeoff Weight (W0) = ', num2str(W0), ' lb']);

% Now you can compute actual empty weight:
We = We_fraction * W0;
disp(['Estimated Empty Weight (We) = ', num2str(We), ' lb']);

% Fuel weight:
Wf = fuel_fraction * W0;
disp(['Estimated Fuel Weight (Wf) = ', num2str(Wf), ' lb']);

disp('---------------------');