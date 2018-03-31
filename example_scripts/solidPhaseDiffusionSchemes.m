% LIONSIMBA example script
% solidPhaseDiffusionSchemes: this script presents a comparison of the simulation results between
% a spectral scheme and a finite difference method when used for the resolution of the solid phase diffusion equations

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
% Define the parameters structure.
param{1} = Parameters_init;

% Use Fick's law of diffusion for solid particles
param{1}.SolidPhaseDiffusion=3;

% Use Finite Difference scheme for the resolution of the solid diffusion equation
param{1}.SolidPhaseDiffusionNumericalScheme =1;

out1 = startSimulation(t0,tf,[],-15,param);

% Use Spectral scheme for the resolution of the solid diffusion equation
param{1}.SolidPhaseDiffusionNumericalScheme =2;

out2 = startSimulation(t0,tf,[],-15,param);

%% Plot the results

figure(1)
plot(out1.time{1},out1.Voltage{1},'LineWidth',3)
hold on
plot(out2.time{1},out2.Voltage{1},'--','LineWidth',3)
xlabel('Time [s]')
ylabel('Voltage [V]')
legend('FDM','SM')
grid on

figure(2)
plot(out1.time{1},out1.Temperature{1},'LineWidth',3)
hold on
plot(out2.time{1},out2.Temperature{1},'--','LineWidth',3)
xlabel('Time [s]')
ylabel('Temperature [K]')
legend('FDM','SM')
grid on

figure(3)
plot(out1.time{1},out1.cs_surface{1}(:,[1,10,15,20]),'LineWidth',3)
hold on
plot(out2.time{1},out2.cs_surface{1}(:,[1,10,15,20]),'--','LineWidth',3)
xlabel('Time [s]')
ylabel('Solid phase concentration [mol/m^3]')
legend('FDM','SM')
grid on