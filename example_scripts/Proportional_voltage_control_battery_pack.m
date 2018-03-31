% LIONSIMBA example script
% Proportional_voltage_control_battery_pack: simulates the usage of a proportional
% controller in order to carry out a charging of a battery composed of 2 series-connected
% cells

%   This file is part of the LIONSIMBA Toolbox
%
%	Official web-site: 	http://sisdin.unipv.it/labsisdin/lionsimba.php
% 	Official GitHUB: 	https://github.com/lionsimbatoolbox/LIONSIMBA
%
%   LIONSIMBA: A Matlab framework based on a finite volume model suitable for Li-ion battery design, simulation, and control
%   Copyright (C) 2016-2018 :Marcello Torchio, Lalo Magni, Davide Raimondo,
%                            University of Pavia, 27100, Pavia, Italy
%                            Bhushan Gopaluni, Univ. of British Columbia, 
%                            Vancouver, BC V6T 1Z3, Canada
%                            Richard D. Braatz, 
%                            Massachusetts Institute of Technology, 
%                            Cambridge, Massachusetts 02142, USA
%   
%   Main code contributors to LIONSIMBA 2.0:
%                           Ian Campbell, Krishnakumar Gopalakrishnan,
%                           Imperial college London, London, UK
%
%   LIONSIMBA is a free Matlab-based software distributed with an MIT
%   license.

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
param{1}.OperatingMode = 4;

% Define the script which will provide the applied current density values.
% Check the file mentioned here below to understand how the proportional
% action has been defined
param{1}.CurrentDensityFunction = @getPcontrolCurrentPack;

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
plot(results.time{1},results.curr_density,'LineWidth',6)
box on
grid on
xlim([0 results.time{1}(end)])
ylim([min(results.curr_density) max(results.curr_density)])
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