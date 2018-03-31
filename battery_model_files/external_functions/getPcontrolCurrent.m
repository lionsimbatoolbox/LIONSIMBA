function I = getPcontrolCurrent(t,t0,tf,x,param,extra)
%	getPcontrolCurrent returns the value of the input current as a function of
%	the time, states and parameters
%
%       I = getPcontrolCurrent(t,t0,tf,y,param,extra)
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

% Get the value of the Voltage out of the current battery states.
V = x(param{1}.Phis_indices(1))-x(param{1}.Phis_indices(end));

% Please note the usage of fields like param{1}.Phis_indices etc. These are
% arrays evaluated during the call to startSimulation and they store the
% relative position of the variables inside the overall array x. The script
% computeVariablesIndices.m in the simulator tools folder explains how such
% indices are evaluated, and the name of the associated variables. The
% example getPcontrolCurrentPack.m explains how to deal when multiple
% cells are present.

% Define the proportional action. Do not put extreme values, the simulator
% could crash
Kp = 100;

% Define the Voltage Setpoint
V_ref = 4.15;
% Define your linear or nonlinear function of t for evaluate the value
% of I.
I = Kp*(V_ref-V);

end