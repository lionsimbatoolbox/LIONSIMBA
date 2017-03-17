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
% Custom current profile: this script shows how to run a custom current
% profile charge/discharge

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

param{1}.hcell = 1;
param{1}.Tref = 298.15;

param{1}.AbsTol        = 1e-6;
param{1}.RelTol        = 1e-6;

param{1}.CutoffSOC     = 20;

param{1}.SolidPhaseDiffusion = 1;

if(param{1}.SolidPhaseDiffusion == 3)
    multp = param{1}.Nr_p;
    multn = param{1}.Nr_n;
else
    multp = 1;
    multn = 1;
end
I1C = 29.5;

C_rate = 1.5;
%% Discharge section
% Discharge the battery to the 20% of SOC
out = startSimulation(t0,tf,initialState,-25,param);

% Store the Jacobian matrix
param{1}.JacobianFunction = out.JacobianFun;

% Update the initial states
initialState = out.initialState;

% Redefine the cutoff SOC in order to avoid simulation interruptions
param{1}.CutoffSOC     = 2;

% Rest the battery with no current applied
out2 = startSimulation(0,5000,initialState,0,param);

% Update the initial states
initialState = out2.initialState;

% Apply a custom current profile
param{1}.AppliedCurrent = 2;

% Run the simulation
out3 = startSimulation(0,5000,initialState,0,param);

% Update the initial states
initialState = out3.initialState;

% Go back to galvanostatic operating conditions
param{1}.AppliedCurrent = 1;

% Run the simulation
out4 = startSimulation(0,5000,initialState,0,param);
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
plot(time,[out.appliedCurrent;out2.appliedCurrent;out3.appliedCurrent;out4.appliedCurrent],'LineWidth',6)
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