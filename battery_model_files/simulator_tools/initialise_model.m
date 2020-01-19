function [init_point, n_alg, initial_terminal_voltage] = initialise_model(param)
% initialise_model performs analytical initialisation of all variables of interest

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

sign_input_density = evaluate_sign_input_density(param); % evaluates sign of input current/power density as per operating mode

T = param.T_init * ones(param.Nsum+param.Nal+param.Ncu,1);

cs_star_p   = param.cs_p_init*ones(param.Np, 1);
cs_star_n   = param.cs_n_init*ones(param.Nn, 1);
cs_star     = [cs_star_p; cs_star_n];

[Phis_p_init,~,Phis_n_init,~] = param.OpenCircuitPotentialFunction(cs_star,T,param,sign_input_density);

Phis_init = [Phis_p_init; Phis_n_init];

if param.OperatingMode==1 || param.OperatingMode==4
    I_density           = param.I_density;
elseif param.OperatingMode==2 || param.OperatingMode==5
    I_density           = param.P_density/(Phis_init(1)-Phis_init(end)); % No need for linear interpolation since starting from equilibrium
else
    I_density = 0; % dummy value for CV mode
end

ce_init = param.ce_init*[ones(param.Np,1);ones(param.Ns,1);ones(param.Nn,1)];
Phie_init = zeros(param.Np + param.Ns + param.Nn, 1); % designated as the ground potential for this system

solverFlux  = zeros(param.Np+param.Nn, 1);
film = zeros(param.Nn, 1);

jflux_init = ionicFlux(ce_init, cs_star, Phis_init, Phie_init, T, solverFlux, film, param,sign_input_density,I_density);

if param.EnableAgeing == 1
    js_init = 0.483e-5*ones(param.Nn, 1); % totally arbitrary/random value for side-reaction flux
else
    js_init = zeros(param.Nn, 1);
end

init_point = [
    jflux_init;...
    Phis_init;...
    Phie_init;...
    js_init;...
    I_density
    ];

n_alg = length(init_point);
initial_terminal_voltage = diff([Phis_n_init(end);Phis_p_init(1)]); % No need for linear interpolation since starting from equilibrium
end
