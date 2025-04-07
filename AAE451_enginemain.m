close all
clear

T_TO_mil = 76.3; % mil standard takeoff thrust (kN)
T_TO_AB = 131.2; % full afterburner takeoff thrust (kN)
correct_mil = 1; % SFC correction factor for dry thrust
correct_AB = 1; % SFC correction factor for wet thrust
iter = 1;

M = 0.28:0.05:0.9;
h = 0:5000:35000;
% M = 0.9;
% h = [0 15000 35000];

[SFC_AB, T_vel, T_mil, T_AB, correct_mil, correct_AB, iter] = AAE451_Thrust(0.01,0,T_TO_mil,T_TO_AB,correct_mil,correct_AB,iter);

for i = 1:length(M)
    for j = 1:length(h)
    [SFC_AB(i,j), T_vel(i,j), T_mil(i,j), T_AB(i,j), correct_mil, correct_AB, iter] = AAE451_Thrust(M(i),h(j),T_TO_mil,T_TO_AB,correct_mil, correct_AB,iter);
    end
end

% clean up T_AB data
temp = find(isfinite(T_AB(:,1)));

figure(1)
plot(M,T_AB)
for j = 1:length(h)
    if(j==1)
        plot([0 M(temp)],[29500; T_AB(temp,j)],'LineWidth',1.5)
        hold on
        txt = ['h = ',num2str(h(j)), ' ft'];
        text(M(temp(1)),T_AB(temp(1),j),txt,'HorizontalAlignment','left','VerticalAlignment','top')
    end
    plot(M(temp),T_AB(temp,j),'LineWidth',1.5)
    txt = ['h = ',num2str(h(j)), ' ft'];
    text(M(temp(1)),T_AB(temp(1),j),txt,'HorizontalAlignment','left','VerticalAlignment','top')
end
grid on
xlabel('Mach Number')
ylabel('Thrust (lbf)')
title('Thrust Variation With Full Afterburner')
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.TickLabelFormat = '%,g';
ylim([0 30000])
xlim([0 0.9])

figure(2)
for j = 1:length(h)
    if(j==1)
        plot([0 M],[17155; T_mil(:,j)],'LineWidth',1.5)
        hold on
        txt = ['h = ',num2str(h(j)), ' ft'];
        text(M(temp(1)),T_mil(temp(1),j),txt,'HorizontalAlignment','left','VerticalAlignment','top')
    end
    plot(M,T_mil(:,j),'LineWidth',1.5)
    txt = ['h = ',num2str(h(j)), ' ft'];
    text(M(temp(1)),T_mil(temp(1),j),txt,'HorizontalAlignment','left','VerticalAlignment','top')
end
grid on
xlabel('Mach Number')
ylabel('Thrust (lbf)')
title('Thrust Variation With Mil Power')
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.TickLabelFormat = '%,g';
ylim([0 17200])
xlim([0 0.9])
