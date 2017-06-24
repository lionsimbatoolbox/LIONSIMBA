%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getPcontrolCurrent returns the value of the input current as a function of
% the time, states and parameters
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

function I = getPcontrolCurrent(t,t0,tf,x,param,extra)

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
V_ref = 4.17;
% Define your linear or nonlinear function of t for evaluate the value
% of I.
I = Kp*(V_ref-V);

end