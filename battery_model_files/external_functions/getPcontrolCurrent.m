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
%               - y         : contains the array of all the states
%                             (differential and algebraic) at time t
%               - param     : contains the parameters structure
%               - extra     : extra parameters
%       Outputs:
%               - I     : Applied current desnity [A/m^2]

function I = getPcontrolCurrent(t,t0,tf,y,param,extra)

% Get the value of the Voltage out of the current battery states
V = y(param.Phis_indices(1))-y(param.Phis_indices(end));

% Define the proportional action. Do not put extreme values, the simulator
% could crash
Kp = 100;

% Define the Voltage Setpoint
V_ref = 4.17;
% Define your linear or nonlinear function of t for evaluate the value
% of I.
I = Kp*(V_ref-V);

end