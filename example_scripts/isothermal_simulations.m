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
% Isothermal scenario: this script provides the example simulation shown in the
% paper.

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
% LIONSAMBA will stop automatically when reached the Cutoff Voltage of
% 2.5V.

out1 = startSimulation(t0,tf,[],-15,param);

% Store the Jacobian matrix
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