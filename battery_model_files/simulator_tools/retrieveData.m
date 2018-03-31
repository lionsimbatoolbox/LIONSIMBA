function [ce_t, cs_barrato_t, T_t, jflux_t, Phis_t, Phie_t, cs_star_t, SOC_t, film_t, js_t, Q_t,t_tot] = retrieveData(ce_t, cs_barrato_t, T_t, jflux_t, Phis_t, Phie_t,cs_star_t, SOC_t, film_t, js_t, Q_t,t_tot, y, t, param)
%	retrieveData retrieves the data which is returned in the results structure.

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

% Extract differential variables after the integration process
ce_t            = [ce_t;y(param.ce_indices)];
cs_barrato_t    = [cs_barrato_t;y(param.cs_average_indices)];
T_t             = [T_t;y(param.T_indices)];
film_t          = [film_t;y(param.film_indices)];
Q_t             = [Q_t;y(param.Q_indices)];
% Extract the algebraic variables after the integration process
jflux_t         = [jflux_t;y(param.jflux_indices)];
Phis_t          = [Phis_t;y(param.Phis_indices)];
Phie_t          = [Phie_t;y(param.Phie_indices)];
js_t            = [js_t;y(param.js_indices)];

% Check if Fick's law of diffusion is used. This is required to define the
% correct way how to evaluate the SOC.
if(param.SolidPhaseDiffusion~=3)
    cs_average = sum(cs_barrato_t(end,param.Np+1:end))/(param.Nn);  % cs_average in neg electrode
else
    cs_average = sum(cs_barrato_t(end, (param.Np*param.Nr_p) +1:end))/(param.Nn*param.Nr_n); % cs_average in neg electrode
end

Sout = 100*((cs_average/param.cs_maxn) - param.theta_min_neg)/(param.theta_max_neg- param.theta_min_neg); % cell-soc percent

cs_star_t       = [cs_star_t;surfaceConcentration(cs_barrato_t(end,:)',jflux_t(end,:)',Q_t(end,:)',T_t(end,:)',param)'];
SOC_t           = [SOC_t;Sout];
t_tot           = [t_tot;t];
end
