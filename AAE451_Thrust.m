function [SFC_mil, T_vel, T_mil, T_AB,correct_mil, correct_AB, iter] = AAE451_Thrust(M,h,T_TO_mil,T_TO_AB,correct_mil, correct_AB, iter)

SFC_mil = AAE451_SFCmodel(M,h,T_TO_mil);
SFC_AB = AAE451_SFCmodel(M,h,T_TO_AB);

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

%% calculate thrust at desired altitude and mach number
ft_to_m = 0.3048; % conversion from feet to meters for standard atmosphere calculator
h_metric = h*ft_to_m;

ms_to_knots = 1.9438444924; % converts m/s to knots (or KTAS)

[T, a, P, rho] = atmosisa(h_metric);

a = a*ms_to_knots; % speed of sound (KTAS)
velocity = M*a; % cruise velocity (KTAS)

% interpolate at velocity
temp_greater = find(mdot_20(:,1)>velocity);
temp_less = find(mdot_20(:,1)<velocity);

% interpolate at altitude
mdot_h = zeros(1,length(mdot_20(:,2)));
if(h>=0 && h<= 4000)
    for i = 1:length(mdot_20(:,2))
        mdot_h(i) = interp1([0 4000], ...
        [mdot_20(i,2) mdot_20(i,3)],h,'linear');
    end
elseif(h>4000 && h<=8000)
    for i = 1:length(mdot_20(:,2))
        mdot_h(i) = interp1([4000 8000], ...
        [mdot_20(i,3) mdot_20(i,4)],h,'linear');
    end
elseif(h>8000 && h<=12000)
    for i = 1:length(mdot_20(:,2))
        mdot_h(i) = interp1([8000 12000], ...
        [mdot_20(i,4) mdot_20(i,5)],h,'linear');
    end
elseif(h>12000 && h<=16000)
    for i = 1:length(mdot_20(:,2))
        mdot_h(i) = interp1([12000 16000], ...
        [mdot_20(i,5) mdot_20(i,6)],h,'linear');
    end
elseif(h>16000 && h<=20000)
    for i = 1:length(mdot_20(:,2))
        mdot_h(i) = interp1([16000 20000], ...
        [mdot_20(i,6) mdot_20(i,7)],h,'linear');
    end
    elseif(h>20000 && h<=25000)
    for i = 1:length(mdot_20(:,2))
        mdot_h(i) = interp1([20000 25000], ...
        [mdot_20(i,7) mdot_20(i,8)],h,'linear');
        if(isnan(mdot_h(i)))
            mdot_h(i) = mdot_20(i,7);
        end
    end
    elseif(h>25000 && h<=30000)
    for i = 1:length(mdot_20(:,2))
        mdot_h(i) = interp1([25000 30000], ...
        [mdot_20(i,8) mdot_20(i,9)],h,'linear');
        if(isnan(mdot_h(i)))
            mdot_h(i) = mdot_20(i,8);
        end
    end
    elseif(h>30000 && h<=35000)
    for i = 1:length(mdot_20(:,2))
        mdot_h(i) = interp1([30000 35000], ...
        [mdot_20(i,9) mdot_20(i,10)],h,'linear');
        if(isnan(mdot_h(i)))
            mdot_h(i) = mdot_20(i,9);
        end
    end
end
mdot_h = transpose(mdot_h);

if(M>0.28)
mdot_vel = interp1([mdot_20(temp_less(end),1) mdot_20(temp_greater(1),1)], ...
     [mdot_h(temp_less(end)) mdot_h(temp_greater(1))],velocity,'linear');

T_vel = mdot_vel/(SFC_mil*correct_mil);
else
    T_vel = NaN;
end

T_mil = mdot_h(end)/(SFC_mil*correct_mil);
if(M==0.01 && h==0)
    correct_mil = T_mil/17000;
end

clear mdot_h
%% afterburner

if(M<0.4)
    M = 0.45;
end
% interpolate on altitude curves about mach number
if(h>=0 && h<= 5000)
    temp_greater1 = find(mdot_AB_SL(:,1)>M);
    temp_less1 = find(mdot_AB_SL(:,1)<M);
    temp_greater2 = find(mdot_AB_5(:,1)>M);
    temp_less2 = find(mdot_AB_5(:,1)<M);

    mdot_h1 = interp1([mdot_AB_SL(temp_less1(end),1) mdot_AB_SL(temp_greater1(1),1)],...
        [mdot_AB_SL(temp_less1(end),2) mdot_AB_SL(temp_greater1(1),2)],M);
    mdot_h2 = interp1([mdot_AB_5(temp_less2(end),1) mdot_AB_5(temp_greater2(1),1)],...
        [mdot_AB_5(temp_less2(end),2) mdot_AB_5(temp_greater2(1),2)],M);
    weight = h/5000;
    mdot_h = mdot_h1*(1-weight) + mdot_h2*weight;
elseif(h>5000 && h<=10000)
    temp_greater1 = find(mdot_AB_5(:,1)>M);
    temp_less1 = find(mdot_AB_5(:,1)<M);
    temp_greater2 = find(mdot_AB_10(:,1)>M);
    temp_less2 = find(mdot_AB_10(:,1)<M);

    mdot_h1 = interp1([mdot_AB_5(temp_less1(end),1) mdot_AB_5(temp_greater1(1),1)],...
        [mdot_AB_5(temp_less1(end),2) mdot_AB_5(temp_greater1(1),2)],M);
    mdot_h2 = interp1([mdot_AB_10(temp_less2(end),1) mdot_AB_10(temp_greater2(1),1)],...
        [mdot_AB_10(temp_less2(end),2) mdot_AB_10(temp_greater2(1),2)],M);
    weight = h/10000;
    mdot_h = mdot_h1*(1-weight) + mdot_h2*weight;
elseif(h>10000 && h<=15000)
    temp_greater1 = find(mdot_AB_10(:,1)>M);
    temp_less1 = find(mdot_AB_10(:,1)<M);
    temp_greater2 = find(mdot_AB_15(:,1)>M);
    temp_less2 = find(mdot_AB_15(:,1)<M);

    mdot_h1 = interp1([mdot_AB_10(temp_less1(end),1) mdot_AB_10(temp_greater1(1),1)],...
        [mdot_AB_10(temp_less1(end),2) mdot_AB_10(temp_greater1(1),2)],M);
    mdot_h2 = interp1([mdot_AB_15(temp_less2(end),1) mdot_AB_15(temp_greater2(1),1)],...
        [mdot_AB_15(temp_less2(end),2) mdot_AB_15(temp_greater2(1),2)],M);
    weight = h/15000;
    mdot_h = mdot_h1*(1-weight) + mdot_h2*weight;
elseif(h>15000 && h<=20000)
    temp_greater1 = find(mdot_AB_15(:,1)>M);
    temp_less1 = find(mdot_AB_15(:,1)<M);
    temp_greater2 = find(mdot_AB_20(:,1)>M);
    temp_less2 = find(mdot_AB_20(:,1)<M);

    mdot_h1 = interp1([mdot_AB_15(temp_less1(end),1) mdot_AB_15(temp_greater1(1),1)],...
        [mdot_AB_15(temp_less1(end),2) mdot_AB_15(temp_greater1(1),2)],M);
    mdot_h2 = interp1([mdot_AB_20(temp_less2(end),1) mdot_AB_20(temp_greater2(1),1)],...
        [mdot_AB_20(temp_less2(end),2) mdot_AB_20(temp_greater2(1),2)],M);
    weight = h/20000;
    mdot_h = mdot_h1*(1-weight) + mdot_h2*weight;
elseif(h>20000 && h<=25000)
    temp_greater1 = find(mdot_AB_20(:,1)>M);
    temp_less1 = find(mdot_AB_20(:,1)<M);
    temp_greater2 = find(mdot_AB_25(:,1)>M);
    temp_less2 = find(mdot_AB_25(:,1)<M);

    mdot_h1 = interp1([mdot_AB_20(temp_less1(end),1) mdot_AB_20(temp_greater1(1),1)],...
        [mdot_AB_20(temp_less1(end),2) mdot_AB_20(temp_greater1(1),2)],M);
    mdot_h2 = interp1([mdot_AB_25(temp_less2(end),1) mdot_AB_25(temp_greater2(1),1)],...
        [mdot_AB_25(temp_less2(end),2) mdot_AB_25(temp_greater2(1),2)],M);
    weight = h/25000;
    mdot_h = mdot_h1*(1-weight) + mdot_h2*weight;
elseif(h>25000 && h<=30000)
    temp_greater1 = find(mdot_AB_25(:,1)>M);
    temp_less1 = find(mdot_AB_25(:,1)<M);
    temp_greater2 = find(mdot_AB_30(:,1)>M);
    temp_less2 = find(mdot_AB_30(:,1)<M);

    mdot_h1 = interp1([mdot_AB_25(temp_less1(end),1) mdot_AB_25(temp_greater1(1),1)],...
        [mdot_AB_25(temp_less1(end),2) mdot_AB_25(temp_greater1(1),2)],M);
    mdot_h2 = interp1([mdot_AB_30(temp_less2(end),1) mdot_AB_30(temp_greater2(1),1)],...
        [mdot_AB_30(temp_less2(end),2) mdot_AB_30(temp_greater2(1),2)],M);
    weight = h/30000;
    mdot_h = mdot_h1*(1-weight) + mdot_h2*weight;
elseif(h>30000 && h<=35000)
    temp_greater1 = find(mdot_AB_30(:,1)>M);
    temp_less1 = find(mdot_AB_30(:,1)<M);
    temp_greater2 = find(mdot_AB_35(:,1)>M);
    temp_less2 = find(mdot_AB_35(:,1)<M);

    mdot_h1 = interp1([mdot_AB_30(temp_less1(end),1) mdot_AB_30(temp_greater1(1),1)],...
        [mdot_AB_30(temp_less1(end),2) mdot_AB_30(temp_greater1(1),2)],M);
    mdot_h2 = interp1([mdot_AB_35(temp_less2(end),1) mdot_AB_35(temp_greater2(1),1)],...
        [mdot_AB_35(temp_less2(end),2) mdot_AB_35(temp_greater2(1),2)],M);
    weight = h/35000;
    mdot_h = mdot_h1*(1-weight) + mdot_h2*weight;
end


T_AB = mdot_h*1000/(SFC_AB*correct_AB);
if(M==0.45 && h==0 && iter == 1)
    correct_AB = T_AB/37000;
end
if(M==0.45)
    T_AB = NaN;
end

iter = iter + 1;

end
