%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function is used to evaluate the Jacobian Matrix of the P2D model
% according to the specifications of the IDA numerical solver. Please refer
% to the IDS user's guide for additional information on this function
%
%
function [J, flag, new_data] = jacobianFunction(t, x, xp, rr, cj, data)

% Extract the function object previously obtained using CasADi
fJ          = data.fJ;

% Evaluate the Jacobian with respect to the current values of the states
% and their time derivatives.
J           = full(fJ(x,cj));

% Return the values
flag        = 0;
new_data    = [];
end

