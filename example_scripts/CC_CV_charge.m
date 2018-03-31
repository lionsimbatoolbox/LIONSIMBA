% LIONSIMBA example script
% CC_CV_charge scenario: this script provides an example to run a
% constant current-constant voltage charge of a cell.

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
close all

%% Parameters
% Define the integration times.
t0 = 0;
tf = 10^4;

% Define the initial state structure
initialState.Y  = [];
initialState.YP = [];

% Define the parameters structure.
param{1}               = Parameters_init;

param{1}.Np            = 10;
param{1}.Ns            = 10;
param{1}.Nn            = 10;

param{1}.hcell         = 1;
param{1}.Tref          = 298.15;

param{1}.AbsTol        = 1e-6;
param{1}.RelTol        = 1e-6;

param{1}.CutoffSOC     = 20;

param{1}.SolidPhaseDiffusion = 3;

if(param{1}.SolidPhaseDiffusion == 3)
    multp = param{1}.Nr_p;
    multn = param{1}.Nr_n;
else
    multp = 1;
    multn = 1;
end

% Define 1C rate of current
I1C = 29.23;

% Run with a 1.5C-rate in the CC stage
C_rate = 1.5;
%% Discharge section
% Discharge the battery to the 20% of SOC
out = startSimulation(t0,tf,initialState,-25,param);

% Store the Jacobian Matrix
param{1}.JacobianFunction = out.JacobianFun;

initialState = out.initialState;

% Redefine the cutoff SOC in order to avoid simulation interruptions
param{1}.CutoffSOC     = 2;

% Rest the battery with no current applied
out2 = startSimulation(0,5000,initialState,0,param);

initialState = out2.initialState;
%% Constant Current
% Set the cutover Voltage to 4.17 V in order to trigger the event for CC-CV
% switching
param{1}.CutoverVoltage = 4.17;

% Charge the battery at 1.5C with CC protocol and stops when CutoverVoltage
% is reached
out3 = startSimulation(t0,tf,initialState,C_rate*I1C,param);

initialState = out3.initialState;
%% Constant Voltage
% Define a new parameter strutcutre for CV cycling
param2 = param;

% Empty the Jacobian matrix!! During the CV stage, the set of equations
% changes
param2{1}.JacobianFunction = [];

% Potentiostatic operations
param2{1}.OperatingMode = 3;
% Change the cutover voltage to avoid simulation interruptions
param2{1}.CutoverVoltage = 4.3;
% Potentiostatic reference
param2{1}.V_reference = out3.Voltage{1}(end);

% Run the simulation in CV
out4 = startSimulation(t0,10000,initialState,0,param2);
initialState = out4.initialState;
%% Plot the results

% Concatenate all the time results
time = [out.time{1};out2.time{1}+out.time{1}(end);out3.time{1}+out2.time{1}(end)+out.time{1}(end);out4.time{1}+out3.time{1}(end)+out2.time{1}(end)+out.time{1}(end)];

figure(1)
plot(time,[out.Voltage{1};out2.Voltage{1};out3.Voltage{1};out4.Voltage{1}],'LineWidth',6)
hold on
xlabel('Time [s]')
ylabel('Voltage [V]')
grid on
box on
title('Cell Voltage')

figure(2)
plot(time,[out.SOC{1};out2.SOC{1};out3.SOC{1};out4.SOC{1}],'LineWidth',6)
hold on
xlabel('Time [s]')
ylabel('SOC [%]')
grid on
box on
title('Cell SOC')

figure(3)
plot(time,[out.Temperature{1}(:,end);out2.Temperature{1}(:,end);out3.Temperature{1}(:,end);out4.Temperature{1}(:,end)],'LineWidth',6)
hold on
xlabel('Time [s]')
ylabel('Temperature[K]')
grid on
box on
title('Cell Temperature')

figure(4)
plot(time,[out.curr_density;out2.curr_density;out3.curr_density;out4.curr_density],'LineWidth',6)
hold on
xlabel('Time [s]')
ylabel('Applied Current [A/m^2]')
grid on
box on
title('Cell Input Current')

figure(5)
plot(time,[out.cs_average{1}(:,[1 param{1}.Np*multp]);out2.cs_average{1}(:,[1 param{1}.Np*multp]);out3.cs_average{1}(:,[1 param{1}.Np*multp]);out4.cs_average{1}(:,[1 param{1}.Np*multp])],'LineWidth',6)
hold on
line([time(1) time(end)],[param{1}.cs_maxp param{1}.cs_maxp],'LineWidth',6,'LineStyle','--')
xlabel('Time [s]')
ylabel('Positive electrode average concentration [mol/m^3]')
grid on
box on
legend('c^* x=0','c^* x=l_p','c^* max')

figure(6)
plot(time,[out.cs_average{1}(:,[(param{1}.Np+1)*multn (param{1}.Np+param{1}.Nn)*multn]);out2.cs_average{1}(:,[param{1}.Np+1 param{1}.Np+param{1}.Nn]*multn);out3.cs_average{1}(:,[param{1}.Np+1 param{1}.Np+param{1}.Nn]*multn);out4.cs_average{1}(:,[param{1}.Np+1 param{1}.Np+param{1}.Nn]*multn)],'LineWidth',6)
hold on
line([time(1) time(end)],[param{1}.cs_maxn param{1}.cs_maxn],'LineWidth',6,'LineStyle','--')
xlabel('Time [s]')
ylabel('Negative electrode average concentration [mol/m^3]')
grid on
box on
legend('c^* x=l_p+l_s','c^* x=l_p+l_s+l_n','c^* max')
