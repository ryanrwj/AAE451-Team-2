%% HDI Drag Polar Script
%  Author: Team 2
%  Date:   3/2/25

W0 = 21184;             % lb, MTOW
W0_S = 92;              % lb/ft^2, wing loading
S = W0/W0_S;            % ft^2, wing area
S_wet = S*2*1.02;       % ft^2, wetted area of wing
AR = 3.2;               % wing aspect ratio
b = sqrt(AR*S_wet);     % ft, wingspan
n_missiles = 0:2:8;     % number of missiles
C_D0missile = 0.003;    % drag added by missile
e0 = 0.8;               % Oswald efficiency factor

C_fe = 0.0035;          % equivalent skin friction coefficient (Raymer 12.3)
C_D0clean = C_fe*S_wet/S;  
K = 1/(pi*AR*e0);

% C_D0 = sum(C_fc*FF_c*Q_c*Swet_c)/S + C_Dmisc + C_DLP

C_L = linspace(-0.5, 2);
C_D = zeros(length(C_L), length(n_missiles));

figure; hold on;
for i = 1:5
    C_D0 = C_D0clean + C_D0missile*n_missiles(i);
    C_D(:, i) = C_D0 + K*C_L.^2;
    plot(C_L, C_D(:, i)), grid;
end

legend('0 missiles', '2 missiles', '4 missiles', '6 missiles', '8 missiles')
title('C_D vs. C_L')
xlabel('C_L')
ylabel('C_D')