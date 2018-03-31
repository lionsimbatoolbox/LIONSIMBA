function [x0_alg, n_alg] =   algebraicInitialConditions(param)
% algebraicInitialConditions returns the number of algebraic variables
% DEPRECATED in favour of initialise_model.m

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

% Initial guess for the algebraic variables
jflux_init          = [-0.43e-5*ones(param.Np,1);0.483e-5*ones(param.Nn,1)];
Phis_init           = [4.2*ones(param.Np,1);0.074*ones(param.Nn,1)];
Phie_init           = zeros(param.Nsum,1);
js_init             = 0.483e-5*ones(param.Nn,1);

% Never mind: these variables returned are not used, since we have analytical initialisation!
% DEPRECATED in favor of initialise_model.m (analytical initialisation)
if param.OperatingMode==1 || param.OperatingMode==4
    I_density       = param.I_density;
elseif param.OperatingMode==2 || param.OperatingMode==5
    I_density       = param.P_density/(Phis_init(1)-Phis_init(end)); % no need for linear interpolation since starting from equilibrium
else
    I_density=0; % dummy value for CV mode
end

% Build the array of algebraic initial conditions
x0_alg              = [
    jflux_init;...
    Phis_init;...
    Phie_init;...
    js_init;...
    I_density
    ];

n_alg = length(x0_alg);
end
