function I = getPcontrolCurrentPack(t,t0,tf,x,param,extra)
%	getPcontrolCurrentPack returns the value of the input current as a function of
%	the time, states and parameters
%
%       I = getPcontrolCurrentPack(t,t0,tf,y,param,extra)
%
%       Inputs:
%               - t         : value of the current time step
%               - t0        : initial integration time
%               - tf        : final integration time
%               - x         : contains the array of all the states
%                             (differential and algebraic) at time t
%               - param     : contains the parameters structure
%               - extra     : extra parameters
%       Outputs:
%               - I     : Applied current desnity [A/m^2]

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

% This script provides the value of the applied current density as a
% function of the voltage across the battery pack.

% Firstly we extract the variables related to the first and second cell. To
% this aim the field x_indices of the param structure contains the values
% of the absolute positionong of the variables in the overall x array.
% Indeed, when running simulations with multiple cells, the states array
% (x) will contain as many rows - for a given time instant t - as the sum
% of the states of all the cells involved in the simulation. x_indices
% stores, for each cell, what are the indices of the overall x in which the
% variables of a given cell are stored.

% Extract the first cell variables
cell1_variables = x(param{1}.x_index);
% Extract the second cell variables
cell2_variables = x(param{2}.x_index);

% At this stage, since we have put the variables of the 2 cells into
% cell1_variables and cell2_variables, use their relative indices to
% extract the exact values.

% Cell 1 voltage
V_1 = cell1_variables(param{1}.Phis_indices(1))-cell1_variables(param{1}.Phis_indices(end));

% Cell 2 voltage
V_2 = cell2_variables(param{2}.Phis_indices(1))-cell2_variables(param{2}.Phis_indices(end));

% Define the proportional action. Do not put extreme values, the simulator
% could crash
Kp = 100;

% Define the Voltage Setpoint
V_ref = 8.2;
% Define your linear or nonlinear function of t for evaluate the value
% of I.
I = Kp*(V_ref-(V_1+V_2));

end