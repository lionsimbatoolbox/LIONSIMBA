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


% Clear the workspace
clear all

% Define the parameters structure.
param = Parameters_init;
% Change the CVs number in each battery section

% Positive Current Collector
param.Nal   = 10;
% Positive Electrode
param.Np    = 10;
% Separator
param.Ns    = 10;
% Negative Electrode
param.Nn    = 10;
% Negative Current Collector
param.Nco   = 10;

param.AbsTol = 1e-6;
param.RelTol = 1e-6;

% Define the applied current densities.
input_currents 	= -[30;0;20;0];
% Define the duration of each applied current step.
timings 		= [1500;5000;3000;5000]
% Note that length(timings)==length(input_currents) !!

% param.CutoffSOC = 30;
% 
% out = startSimulation(0,10^5,[],[],-20,param);
% param.CutoffSOC = 1;
% out = startSimulation(0,2000,out.Y,out.YP,0,param);
% Initialize the states and their time derivatives
% Y0 	= out.Y;
% YP0 = out.YP;
Y0 	= [];
YP0 = [];
% Set the initial integration time
t0 	= 0;
tf  = 0;
% Initialize the array where results will be stored
ce_tot 			= [];
Phis_tot 		= [];
t_tot 			= [];
Temperature_tot = [];
SOC_tot         = [];

run_n = 1;
for i=1:length(timings)
    % Set the final integration time
    tf 			= tf + timings(i);
    % Start the simulations
    results 	= startSimulation(t0,tf,Y0,YP0,input_currents(i),param);
    % Concatenate the results
    if(run_n==1)
        ce_tot 			= [ce_tot;results.ce];
        SOC_tot         = [SOC_tot;results.SOC];
        Phis_tot 		= [Phis_tot;results.Phis];
        Temperature_tot = [Temperature_tot;results.Temperature];
        t_tot           = [t_tot;results.time];
        run_n = 2;
    else
        ce_tot 			= [ce_tot;results.ce(2:end,:)];
        SOC_tot         = [SOC_tot;results.SOC(2:end,:)];
        Phis_tot 		= [Phis_tot;results.Phis(2:end,:)];
        Temperature_tot = [Temperature_tot;results.Temperature(2:end,:)];
        t_tot           = [t_tot;results.time(2:end)];
    end
    % Update the initial states.
    Y0 	= results.Y;
    YP0 = results.YP;
    % Update the starting integration time instant.
    t0 		= results.time(end);
    
end

%% Voltage plot

plot(t_tot,Phis_tot(:,1)-Phis_tot(:,end),'LineWidth',6)
hold on
box on
grid on
xlim([0 sum(timings)])
xlabel('Time [s]')
ylabel('Cell Voltage [V]')

%% Input profile plots
figure
hold on
t0 = 0;
tf = 0;
prev_current = input_currents(1);
cur_current = input_currents(2);

for i=1:length(timings)
    if(i<length(timings))
        cur_current = input_currents(i+1);
    end
    tf = tf + timings(i);
    % cur_current = input_currents(i);
    line([t0 tf],[input_currents(i) input_currents(i)],'LineWidth',6)
    line([tf tf],[prev_current cur_current],'LineWidth',6)
    t0 = t0+timings(i);
    prev_current = cur_current;
end
box on
grid on
xlim([0 sum(timings)])
ylim([min(input_currents) max(input_currents)])
xlabel('Time [s]')
ylabel('Applied Current Density [A/m^2]')
%% Temperature plot
figure
plot(t_tot,Temperature_tot(:,end),'LineWidth',6)
box on
grid on
xlim([0 sum(timings)])
xlabel('Time [s]')
ylabel('Temperature [K]')

%% SOC plot
figure
plot(t_tot,SOC_tot,'LineWidth',6)
box on
grid on
xlim([0 sum(timings)])
xlabel('Time [s]')
ylabel('SOC [%]')