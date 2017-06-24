%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIONSAMBA example script
% Proportional_voltage_control_battery_pack: simulates the usage of a proportional
% controller in order to carry out a charging of a battery composed of 2 series-connected
% cells

% Clear the workspace
clear

% Define the integration times.
t0 = 0;
tf = 10^4;

% Define the parameters structure for the battery pack
param{1} = Parameters_init;
param{2} = Parameters_init;

param{1}.CutoffSOC = 20;
% Discharge the cell down to 20% SOC
results = startSimulation(t0,tf,[],-20,param);

% Store the Jacobian matrix for future usage
param{1}.JacobianFunction = results.JacobianFun;

% Push down the CutoffSOC value in order to avoid the premature interruption of
% the simulation
param{1}.CutoffSOC = 2;

% Rest the cell
results = startSimulation(t0,tf,results.initialState,0,param);

% Enable the custom current profile
param{1}.AppliedCurrent = 2;

% Define the script which will provide the applied current density values.
% Check the file mentioned here below to understand how the proportional
% action has been defined
param{1}.CurrentFunction = @getPcontrolCurrentPack;

% Remove the stored Jacobian matrix. This is required because now the
% applied current density is function of the states.
param{1}.JacobianFunction = [];

% Run the proportional control simulation
results = startSimulation(t0,5*tf,results.initialState,0,param);

%% Voltage plot

plot(results.time{1},results.Voltage{1}+results.Voltage{2},'LineWidth',6)
hold on
box on
grid on
xlim([0 results.time{1}(end)])
xlabel('Time [s]')
ylabel('Pack Voltage [V]')

%% Input profile plots
figure
hold on
plot(results.time{1},results.appliedCurrent,'LineWidth',6)
box on
grid on
xlim([0 results.time{1}(end)])
ylim([min(results.appliedCurrent) max(results.appliedCurrent)])
xlabel('Time [s]')
ylabel('Applied Current Density [A/m^2]')
%% Temperature plot
figure
subplot(2,1,1)

plot(results.time{1},results.Temperature{1}(:,end),'LineWidth',6)
title('Temperature [K]')
box on
grid on
xlim([0 results.time{1}(end)])
ylabel('Cell 1')
subplot(2,1,2)
plot(results.time{1},results.Temperature{2}(:,end),'LineWidth',6)
box on
grid on
xlim([0 results.time{2}(end)])
xlabel('Time [s]')
ylabel('Cell 2')

%% SOC plot
figure
subplot(2,1,1)
plot(results.time{1},results.SOC{1},'LineWidth',6)
title('SOC [%]')
box on
grid on
xlim([0 results.time{1}(end)])
xlabel('Time [s]')
ylabel('Cell 1')
subplot(2,1,2)
plot(results.time{1},results.SOC{2},'LineWidth',6)
box on
grid on
xlim([0 results.time{1}(end)])
xlabel('Time [s]')
ylabel('Cell 2')