close all
clear

M = 0.3:0.01:0.9;%linspace(0.3,0.9,100); % mach number 
T_TO_mil = 76.3; % mil-standard takeoff thrust (kN)
T_TO_AB = 131.2; % full afterburner takeoff thrust (kN)
T_h = [518.7 504.4 490.2 475.9 461.7 447.5 429.7 411.9 394.1]*(5/9); % temperature at altitude (K)
SFC_mil = AAE451_SFCmodel(M,T_h,T_TO_mil);
SFC_AB = AAE451_SFCmodel(M,T_h,T_TO_AB);

figure(1)
for i = 1:length(T_h)
    plot(M,SFC_mil(i,:),'LineWidth',1.5)
    hold on
end
legend('Sea Level','4000 feet','8000 feet','12000 feet','16000 feet','20000 feet','25000 feet','30000 feet','35000 feet')
grid on
xlabel('Mach number')
ylabel('SFC (lbm/lbf/hr)')
title('TSFC vs Mach Number Using Mil Thrust')

figure(2)
for i = 1:length(T_h)
    plot(M,SFC_AB(i,:),'LineWidth',1.5)
    hold on
end
legend('Sea Level','4000 feet','8000 feet','12000 feet','16000 feet','20000 feet','25000 feet','30000 feet','35000 feet')
grid on
xlabel('Mach number')
ylabel('SFC (lbm/lbf/hr)')
title('TSFC vs Mach Number Using Wet(with afterburner) Thrust')

data_cruise = readmatrix('Engine Data GE F110-GE-129.txt');
data_AB = readmatrix('GE F110-GE-129 Afterburner Fuel Flow Rate Data.csv');

% cruise fuel flow rate for gross weight of 20000 pounds (lb/hr)
mdot_20 = [data_cruise(4:end,1) data_cruise(4:end,2) data_cruise(4:end,7) data_cruise(4:end,12) data_cruise(4:end,17) ...
    data_cruise(4:end,22) data_cruise(4:end,27) data_cruise(4:end,32) data_cruise(4:end,37) data_cruise(4:end,42)]; 

% first column is mach number, second column is fuel flow in thousands of pounds per hour
mdot_AB_SL = [data_AB(3:end,1) data_AB(3:end,2)]; % sea level
mdot_AB_5 = [data_AB(3:end,3) data_AB(3:end,4)]; % 5000 feet altitude
mdot_AB_10 = [data_AB(3:end,5) data_AB(3:end,6)]; % 10000 feet altitude
mdot_AB_15 = [data_AB(3:end,7) data_AB(3:end,8)]; % 15000 feet altitude
mdot_AB_20 = [data_AB(3:end,9) data_AB(3:end,10)]; % 20000 feet altitude
mdot_AB_25 = [data_AB(3:end,11) data_AB(3:end,12)]; % 25000 feet altitude
mdot_AB_30 = [data_AB(3:end,13) data_AB(3:end,14)]; % 30000 feet altitude
mdot_AB_35 = [data_AB(3:end,15) data_AB(3:end,16)]; % 35000 feet altitude

%% calculations at sea level for thrust at velocity, mil thrust and with afterburner at mach 0.9
soundspeed = 661.479; % SL speed of sound (KTAS)
M_cruise = 0.9;
cruise_vel = M_cruise*soundspeed;

% thrust at velocity
temp_greater = find(mdot_20(:,1)>cruise_vel);
temp_less = find(mdot_20(:,1)<cruise_vel);
mdot_calc = interp1([mdot_20(temp_less(end),1) mdot_20(temp_greater(1),1)], ...
    [mdot_20(temp_less(end),2) mdot_20(temp_greater(1),2)],cruise_vel,'linear');

thrust_vel_SL = mdot_calc/SFC_AB(1,M==M_cruise);

thrust_mil_SL = mdot_20(end,2)/SFC_AB(1,M==M_cruise);

% afterburner
temp_greater = find(mdot_AB_SL(:,1)>M_cruise);
temp_less = find(mdot_AB_SL(:,1)<M_cruise);
mdot_calc = interp1([mdot_AB_SL(temp_less(end),1) mdot_AB_SL(temp_greater(1),1)], ...
    [mdot_AB_SL(temp_less(end),2)*1000 mdot_AB_SL(temp_greater(1),2)*1000],M_cruise,'linear');

thrust_AB_SL = mdot_calc/(SFC_AB(1,M==M_cruise)*2.3);

%% calculations at 15000 feet for thrust at velocity, mil thrust and with afterburner at mach 0.9
soundspeed = 624.035; % SL speed of sound (KTAS)
M_cruise = 0.9;
cruise_vel = M_cruise*soundspeed;

% thrust at velocity
temp_greater = find(mdot_20(:,1)>cruise_vel);
temp_less = find(mdot_20(:,1)<cruise_vel);
mdot_calc = interp1([mdot_20(temp_less(end),1) mdot_20(temp_greater(1),1)], ...
    [mdot_20(temp_less(end),6) mdot_20(temp_greater(1),6)],cruise_vel,'linear');

thrust_vel_15 = mdot_calc/SFC_AB(5,M==M_cruise);

thrust_mil_15 = mdot_20(end,6)/SFC_AB(5,M==M_cruise);

% afterburner
temp_greater = find(mdot_AB_15(:,1)>M_cruise);
temp_less = find(mdot_AB_15(:,1)<M_cruise);
mdot_calc = interp1([mdot_AB_15(temp_less(end),1) mdot_AB_15(temp_greater(1),1)], ...
    [mdot_AB_15(temp_less(end),2)*1000 mdot_AB_15(temp_greater(1),2)*1000],M_cruise,'linear');

thrust_AB_15 = mdot_calc/(SFC_AB(5,M==M_cruise)*2.3);

%% calculations at 35000 feet for thrust at velocity, mil thrust and with afterburner at mach 0.9
soundspeed = 576.419; % SL speed of sound (KTAS)
M_cruise = 0.9;
cruise_vel = M_cruise*soundspeed;

% thrust at velocity
temp_greater = find(mdot_20(:,1)>cruise_vel);
temp_less = find(mdot_20(:,1)<cruise_vel);
mdot_calc = interp1([mdot_20(temp_less(end),1) mdot_20(temp_greater(1),1)], ...
    [mdot_20(temp_less(end),10) mdot_20(temp_greater(1),10)],cruise_vel,'linear');

thrust_vel_35 = mdot_calc/SFC_AB(9,M==M_cruise);

thrust_mil_35 = mdot_20(end,10)/SFC_AB(9,M==M_cruise);

% afterburner
temp_greater = find(mdot_AB_15(:,1)>M_cruise);
temp_less = find(mdot_AB_15(:,1)<M_cruise);
mdot_calc = interp1([mdot_AB_15(temp_less(end),1) mdot_AB_15(temp_greater(1),1)], ...
    [mdot_AB_35(temp_less(end),2)*1000 mdot_AB_35(temp_greater(1),2)*1000],M_cruise,'linear');

thrust_AB_35 = mdot_calc/(SFC_AB(9,M==M_cruise)*2.3);