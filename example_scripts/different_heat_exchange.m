%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIONSAMBA example script
% Different heat exchange coefficients: this script provides the example simulation shown in the
% paper.

% Clear the workspace
clear all

% Define the integration times.
t0 = 0;
tf = 10^4;
% Define the parameters structure. By default the hcell parameter is set to
% 1 [W / (m^2 K)]
param{1} = Parameters_init;

% Change the hcell parameter and set it to 0.01 [W / (m^2 K)]
param{1}.hcell = 0.01;

% Start the simulation. Note that the final integration time is 10^4 and
% LIONSAMBA will stop automatically when reached the Cutoff Voltage of
% 2.5V.

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

