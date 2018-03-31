function [ce_flux_p, ce_flux_ps, ce_flux_s, ce_flux_sn, ce_flux_n] = interpolateElectrolyteConcetrationFluxes(ce,param)
%	interpolateElectrolyteConcetrationFluxes interpolates the electrolyte concentration flux at the edges of control volumes using harmonic mean.

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

% Fluxes within the positive electrode
ce_flux_p = (ce(2:param.Np)-ce(1:param.Np-1))/(param.deltax_p*param.len_p);

% Fluxes at the separator-positive interface
ce_flux_ps = (ce(param.Np+1)-ce(param.Np)) / ((param.deltax_p*param.len_p/2+param.deltax_s*param.len_s/2));

% Fluxes within the separator
ce_flux_s = (ce(param.Np+2:param.Np+param.Ns)-ce(param.Np+1:param.Np+param.Ns-1))/(param.deltax_s*param.len_s);

% Fluxes at the separator-negative interface
ce_flux_sn = (ce(param.Np+param.Ns+1)-ce(param.Np+param.Ns)) / ((param.deltax_n*param.len_n/2+param.deltax_s*param.len_s/2));

% Fluxes within the negative electrode
ce_flux_n = (ce(param.Np+param.Ns+2:end)-ce(param.Np+param.Ns+1:end-1))/(param.deltax_n*param.len_n);

end