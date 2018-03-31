function param = solidPhaseDifferentiationMatrices(param)
%   solidPhaseDifferentiationMatrices preallocates the matrices used for numerical
%	derivatives used for the computation of the Phis index.
%	Given that in the proposed code the solid phase conductivity is considered to be constant across the cell
%	length, the discretization matrices can be assembled only once and
%	used later in the code.

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

% A matrix for the positive electrode
c = ones(param.Np-1,1);
d = -2*ones(param.Np,1);

A_p 				= gallery('tridiag',c,d,c);
A_p(1,1) 			= -1;
A_p(end,end-1:end) 	= [1 -1];

% A matrix for the negative electrode
c = ones(param.Nn-1,1);
d = -2*ones(param.Nn,1);

A_n 			= gallery('tridiag',c,d,c);
A_n(1,1) 		= -1;
A_n(end,end) 	= -1;

% Store the matrices in the param structure for future usage for each
% cell in the pack.
if param.OperatingMode==2 || param.OperatingMode==5 % if operating in applied power density mode
    A_n(end,:)=[]; % Non-linear boundary condition shall be applied in algebraicStates.m file (for the complicated power density BC derivation)
end

param.A_p = A_p;
param.A_n = A_n;
end
