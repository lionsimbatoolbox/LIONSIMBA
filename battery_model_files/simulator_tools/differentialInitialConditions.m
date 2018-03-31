function [cs_average_init, ce_init, T_init, film_init, Q_init, n_diff] =   differentialInitialConditions(param)
%	differentialInitialConditions initializes the values of the differential variables

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

% Check the type of model used for solid diffusion
if(param.SolidPhaseDiffusion==1 || param.SolidPhaseDiffusion==2)
    % This initialization is used when reduced models are employed
    cs_average_init     = [param.cs_p_init*ones(param.Np,1);param.cs_n_init*ones(param.Nn,1)];
elseif (param.SolidPhaseDiffusion==3)
    % If the full model is used (Fick's law), then the initial conditions are
    % modified in order to account for the solid phase diffusion
    % equation structure.
    cs_average_init     = [param.cs_p_init*ones(param.Np*param.Nr_p,1);param.cs_n_init*ones(param.Nn*param.Nr_n,1)];
end
% Initial values for the other differential variables.
ce_init             = param.ce_init*[ones(param.Np,1);ones(param.Ns,1);ones(param.Nn,1)];
T_init              = param.T_init * ones(param.Nsum+param.Nal+param.Ncu,1);
film_init           = zeros(param.Nn,1);
Q_init              = zeros(param.Np+param.Nn,1);

% Store the number of differential variables in each cell.
n_diff= sum([length(cs_average_init) length(ce_init) length(T_init) length(film_init) length(Q_init)]);
end
