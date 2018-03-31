% LIONSIMBA example script
% Isothermal scenario: this script provides the example simulation shown in the
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
clear all

% Define the integration times.
t0 = 0;
tf = 10^4;
% Define the parameters structure.
param{1} = Parameters_init;

% Disable the thermal dynamics
param{1}.TemperatureEnabled = 0;

% Start the simulation. Note that the final integration time is 10^4 and
% LIONSIMBA will stop automatically when reached the Cutoff Voltage of
% 2.5V.

out1 = startSimulation(t0,tf,[],-15,param);

% Store the Jacobian matrix for future computations. This is possible
% because the different scenarios share the same model structure, and the
% only quantity which differs is the applied current density.
param{1}.JacobianFunction = out1.JacobianFun;

% Run the simulation
out2 = startSimulation(t0,tf,[],-30,param);

% Run the simulation
out3 = startSimulation(t0,tf,[],-60,param);

%% Plot the results

figure(1)
plot(out1.time{1},out1.Voltage{1},'LineWidth',6)
hold on
plot(out2.time{1},out2.Voltage{1},'--','LineWidth',6)
plot(out3.time{1},out3.Voltage{1},'-.','LineWidth',6)
xlabel('Time [s]')
ylabel('Voltage [V]')
grid on
box on
