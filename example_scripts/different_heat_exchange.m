% LIONSIMBA example script
% Different heat exchange coefficients: this script provides the example simulation shown in the
% paper.

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
% Define the parameters structure. By default the hcell parameter is set to
% 1 [W / (m^2 K)]
param{1} = Parameters_init;

% Change the hcell parameter and set it to 0.01 [W / (m^2 K)]
param{1}.hcell = 0.01;

% Start the simulation. Note that the final integration time is 10^4 and
% LIONSIMBA will stop automatically when reached the Cutoff Voltage of
% 2.5V. Note that no Jacobian among one simulation and the other. This is
% due to the fact that the model structure changes for each of the
% simulation (because hcell changes) and consequently also the Jacobian
% matrix changes and cannot be used for all the scenarios.

out1 = startSimulation(t0,tf,[],-30,param);

% Change the hcell parameter and set it to 1 [W / (m^2 K)]
param{1}.hcell = 1;

% Run the simulation
out2 = startSimulation(t0,tf,[],-30,param);

% Change the hcell parameter and set it to 100 [W / (m^2 K)]
param{1}.hcell = 100;

% Run the simulation
out3 = startSimulation(t0,tf,[],-30,param);

%% Plot the results

figure(1)
plot(out1.time{1},out1.Temperature{1}(:,end),'LineWidth',6)
hold on
plot(out2.time{1},out2.Temperature{1}(:,end),'--','LineWidth',6)
plot(out3.time{1},out3.Temperature{1}(:,end),'-.','LineWidth',6)
xlabel('Time [s]')
ylabel('Temperature [K]')
grid on
box on

figure(2)
plot(out1.time{1},out1.Voltage{1},'LineWidth',6)
hold on
plot(out2.time{1},out2.Voltage{1},'--','LineWidth',6)
plot(out3.time{1},out3.Voltage{1},'-.','LineWidth',6)
xlabel('Time [s]')
ylabel('Voltage [V]')
grid on
box on

