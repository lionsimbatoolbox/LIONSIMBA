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
% Car Cycling: this script provides the example simulation shown in the
% paper. For simplicity the HEV battery pack is considered composed of only
% one cell and its behavior is plotted after the simulation ends.

% Clear the workspace
clear all

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

% Define the applied current densities.
input_currents = [-29.5;14.75;-14.75;-29.5;-58;-29.5;14.75];
% Define the duration of each applied current step.
timings = [50,10,150,200,5,200,10];
% Note that length(timings)==length(input_currents) !!

% Initialize the states and their time derivatives
initialState.Y      = [];
initialState.YP     = [];
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
    results 	= startSimulation(t0,tf,initialState,input_currents(i),param);
    % Concatenate the results
    if(run_n==1)
        ce_tot 						= [ce_tot;results.ce{1}];
        SOC_tot         			= [SOC_tot;results.SOC{1}];
        Phis_tot 					= [Phis_tot;results.Phis{1}];
        Temperature_tot 			= [Temperature_tot;results.Temperature{1}];
        t_tot           			= [t_tot;results.time{1}];
		param{1}.JacobianFunction 	= results.JacobianFun;
        run_n = 2;
    else
        ce_tot 			= [ce_tot;results.ce{1}(2:end,:)];
        SOC_tot         = [SOC_tot;results.SOC{1}(2:end,:)];
        Phis_tot 		= [Phis_tot;results.Phis{1}(2:end,:)];
        Temperature_tot = [Temperature_tot;results.Temperature{1}(2:end,:)];
        t_tot           = [t_tot;results.time{1}(2:end)];
    end
    % Update the initial states.
    initialState = results.initialState;
    % Update the starting integration time instant.
    t0 		= results.time{1}(end);
    
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