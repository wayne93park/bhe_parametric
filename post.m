clear all; close all; clc;

%% Variable 
rb = [0.14/2, 0.15/2, 0.16/2, 0.17/2, 0.18/2. 0.19/2, 0.20/2]'; % borehole radius [m]
rbs = ["0.07", "0.075", "0.08", "0.085", "0.09", "0.095", "0.1" ]';
%% Load CSVs
COP = csvread('COP.csv')';
dp = csvread('pressure_drop.csv')';
n = csvread('efficiency.csv')';
Tb_mean = csvread('mean_Tb.csv')';
t = csvread('time_cxc_0.1.csv');
annulus = [];
centre = [];

for i = 1:size(rb,1)
    filename_a = 'annulus_pipe_out_temp_cxc_' + rbs(i) + '.csv'
    filename_c = 'centre_pipe_in_temp_cxc_' + rbs(i) + '.csv'
    annulus(:,:,i) = csvread(filename_a);
    centre(:,:,i) = csvread(filename_c);
end

T_in = [];
T_out = [];

for j = 1:size(rb,1)
    T_out(j) = annulus(1,end,j);
    T_in(j) = centre(1,end,j);
end

T_in = T_in';
T_out = T_out';

% Plot
figure(1)
plot(rb,COP)
title('Borehole Radius vs. COP')
xlabel('Rb (m)')
ylabel('COP Total')
grid on

figure(2)
plot(rb,dp)
title('Borehole Radius vs. Pressure Drop')
xlabel('Rb (m)')
ylabel('dP (bar)')
grid on

figure(3)
plot(rb, Tb_mean)
title('Borehole Radius vs. Mean Borehole Wall Temperature')
xlabel('Rb (m)')
ylabel('Tb(C)')
grid on

figure(4)
plot(rb, T_in,'g', rb, T_out, 'b')    
title('Borehole Radius vs. Borehole Inlet/Outlet Temp.')
xlabel('Rb (m)')
ylabel('T (C)')
legend('Inlet Temp.','Outlet Temp.')
grid on

figure(5)
plot(t, annulus(1,:,1), 'r', t, centre(1,:,1), '--r', t,...
    annulus(1,:,3), 'b', t, centre(1,:,3),'--b',t, annulus(1,:,7), 'g',...
    t, centre(1,:,7), '--g')    
title('Borehole Radius  vs. Borehole Inlet/Outlet Temp.')
xlabel('t (day)')
ylabel('T (C)')
legend('T_o (Rb = 0.07 m)','T_i (Rb = 0.07 m)', 'T_o (Rb = 0.08 m)','T_i (Rb = 0.08 m)',...
    'T_o (Rb = 0.1 m)','T_i (Rb = 0.1 m)')
grid on

