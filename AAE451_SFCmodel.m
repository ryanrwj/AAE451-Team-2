function [SFC_h] = AAE451_SFCmodel(M,h,T_TO)

ft_to_m = 0.3048; % conversion from feet to meters for standard atmosphere calculator
h_m = h*ft_to_m;

T0 = 288; % given stagnation temperature (K)
delp = 0.02; % deltaP/P, inlet pressure loss
gamma = 1.4; % ratio of specific heats

altitude = [0 4000 8000 12000 16000 20000 25000 30000 35000]; % altitude increments (ft)
T_h = [518.7 504.4 490.2 475.9 461.7 447.5 429.7 411.9 394.1]*(5/9); % temperature at altitude (K)
%% F110-GE-132 Parameters

OAPR = 30.7; % overall pressure ratio
beta = 0.76; % bypass ratio

%% calcualted parameters

eta_i = 1 - (1.3 + 0.25*beta)*delp; % inlet efficiency
eta_c = -2/(2+T_TO) - 0.1171/(0.1171+beta) - M*0.0541+0.9407; % compressor efficiency
eta_t = -3.403/(3.403+T_TO) + 1.048 - M*0.1553; % turbine efficiency
eta_f = -5.978/(5.978+T_TO) - M*0.1479 - 0.1355/(0.1355+beta) + 1.055; % fan efficiency
eta_n = -2.032/(2.032+T_TO) + 1.008 - M*0.009868; % nozzle efficiency
eta_gasgen = 1 - (0.7*(M^2)*(1 - eta_i))/(1 + 0.2*M^2); % gas generator efficiency

theta = 1 + ((gamma-1)/2)*M^2; % isentropic relation stagnation temperature/static temperature
T_TE = (-8000)/T_TO + 1520; % turbine entry temperature in cruise (K)
phi = T_TE./T_h; % temperature ratio between turbine entry temperature and temp at altitude
chi = theta*(OAPR^((gamma-1)/gamma) - 1);

G = zeros(1,length(T_h));
SFC = zeros(1,length(T_h));
for i = 1:length(T_h)
    G(i) = (phi(i) - (chi/eta_c))*(1 - 1.01/ ( (eta_gasgen^((gamma-1)/gamma))*(chi + theta)*...
    (1 - chi/(phi(i)*eta_c*eta_t))));
%% SFC model
    SFC(i) = (0.697*sqrt(T_h(i)/T0)*(phi(i)-theta-chi/eta_c) ) / ...
    ( sqrt(5*eta_n*(1+eta_f*eta_t*beta)*(G(i)+0.2*(M^2)*beta*eta_c/(eta_f*eta_t))) - ...
    M*(1+beta) );
end

T = atmosisa(h_m,extended=true);
temp_greater = find(T_h>T);
temp_less = find(T_h<T);
if(isempty(temp_less))
    SFC_h = SFC(temp_greater(end));
else
    SFC_h = interp1([T_h(temp_less) T_h(temp_greater)], ...
    [SFC(temp_less) SFC(temp_greater)],T,'linear');
end
