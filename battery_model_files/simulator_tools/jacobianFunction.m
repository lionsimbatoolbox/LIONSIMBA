function [J, flag, new_data] = jacobianFunction(~, x, ~, ~, cj, data)
%	jacobianFunction is used to evaluate the Jacobian Matrix of the P2D model
%	according to the specifications of the IDA numerical solver. Please refer
%	to the IDS user's guide for additional information on this function

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

% Extract the function object previously obtained using CasADi
fJ          = data.fJ;

% Evaluate the Jacobian with respect to the current values of the states
% and their time derivatives.
J           = full(fJ(x,cj));

% Return the values
flag        = 0;
new_data    = [];
end