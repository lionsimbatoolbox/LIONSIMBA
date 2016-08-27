%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GETINPUTCURRENT  returns the value of the input current as a function of
% the time.
%
%       I = GETINPUTCURRENT(t,t0,tf,extra)
%
%       Inputs:
%               - t         : value of the current time step
%               - t0        : initial integration time
%               - tf        : final integration time
%               - extra     : extra parameters
%       Outputs:
%               - I     : Applied current desnity [A/m^2]

function I = getInputCurrent(t,t0,tf,extra)

% Define your linear or nonlinear function of t for evaluate the value
% of I.
% I = (t-t0)/(tf-t0) *(-30) + 0;
I = 30*sin(t/100);
end