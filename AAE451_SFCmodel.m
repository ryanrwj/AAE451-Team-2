function [SFC] =  AAE451_SFCmodel(M,T_h, T_TO)

T0 = 288; % given stagnation temperature (K)
delp = 0.02; % deltaP/P, inlet pressure loss
gamma = 1.4; % ratio of specific heats

altitude = [0 4000 8000 12000 16000 20000 25000 30000 35000]; % altitude increments (ft)

%% F110-GE-132 Parameters

OAPR = 30.7; % overall pressure ratio
%T_TO = 76.3; % take-off thrust (mil-standard) (kN)
beta = 0.76; % bypass ratio

%% calcualted parameters

eta_i = 1 - (1.3 + 0.25*beta)*delp; % inlet efficiency
eta_c = -2/(2+T_TO) - 0.1171/(0.1171+beta) - M*0.0541+0.9407; % compressor efficiency
eta_t = -3.403/(3.403+T_TO) + 1.048 - M*0.1553; % turbine efficiency
eta_f = -5.978/(5.978+T_TO) - M*0.1479 - 0.1355/(0.1355+beta) + 1.055; % fan efficiency
eta_n = -2.032/(2.032+T_TO) + 1.008 - M*0.009868; % nozzle efficiency
eta_gasgen = 1 - (0.7*(M.^2)*(1 - eta_i))/(1 + 0.2*M.^2); % gas generator efficiency

theta = 1 + ((gamma-1)/2)*M.^2; % isentropic relation stagnation temperature/static temperature
T_TE = (-8000)/T_TO + 1520; % turbine entry temperature in cruise (K)
phi = T_TE./T_h; % temperature ratio between turbine entry temperature and temp at altitude
chi = theta*(OAPR^((gamma-1)/gamma) - 1);

for i = 1:length(T_h)
    for j = 1:length(M)
        G(i,j) = (phi(i) - (chi(j)./eta_c(j))).*(1 - 1.01./ ( (eta_gasgen^((gamma-1)/gamma)).*(chi(j) + theta(j)).*...
        (1 - chi(j)./(phi(i)*eta_c(j).*eta_t(j)))));
%% SFC model
        SFC(i,j) = (0.697*sqrt(T_h(i)/T0)*(phi(i)-theta(j)-chi(j)/eta_c(j)) ) / ...
        ( sqrt(5*eta_n(j)*(1+eta_f(j)*eta_t(j)*beta)*(G(i,j)+0.2*(M(j)^2)*beta*eta_c(j)/(eta_f(j)*eta_t(j)))) - ...
        M(j)*(1+beta) );

    end
end

end