% LIONSIMBA example script
% Car Cycling: this script provides an alternative way to run 
% the example simulation shown in the paper. For simplicity the HEV battery 
% pack is considered composed of only one cell and its behaviour is plotted after the simulation ends. 
% This is the second of two examples. In this case, only a single
% simulation is run, and the capabilities of the simulator to identify
% discontinuities in the applied input current are used to switch among the
% different values of I.

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

% Define the parameters structure.
param{1} = Parameters_init;
% Change the CVs number in each battery section

% Positive Current Collector
param{1}.Nal   = 10;
% Positive Electrode
param{1}.Np    = 30;
% Separator
param{1}.Ns    = 30;
% Negative Electrode
param{1}.Nn    = 30;
% Negative Current Collector
param{1}.Nco   = 10;


% Note that length(timings)==length(input_currents) !!

% Initialize the states and their time derivatives
initialState.Y      = [];
initialState.YP     = [];
% Set the initial integration time
t0 	= 0;
tf  = 625;
param{1}.OperatingMode     = 4;

% See the getCarCurrent.m file in order to understand how the piecewise
% input currents have to be defined to properly run with LIONSIMBA
param{1}.CurrentDensityFunction    = @getCarCurrent;

% Start the simulations
results 	= startSimulation(t0,tf,initialState,[],param);


%% Voltage plot

plot(results.time{1},results.Phis{1}(:,1)-results.Phis{1}(:,end),'LineWidth',6)
hold on
box on
grid on
xlim([0 results.time{1}(end)])
xlabel('Time [s]')
ylabel('Cell Voltage [V]')

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
plot(results.time{1},results.Temperature{1}(:,end),'LineWidth',6)
box on
grid on
xlim([0 results.time{1}(end)])
xlabel('Time [s]')
ylabel('Temperature [K]')

%% SOC plot
figure
plot(results.time{1},results.SOC{1},'LineWidth',6)
box on
grid on
xlim([0 results.time{1}(end)])
xlabel('Time [s]')
ylabel('SOC [%]')
