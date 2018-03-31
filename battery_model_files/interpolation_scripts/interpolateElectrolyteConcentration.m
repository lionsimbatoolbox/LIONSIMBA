function [ce_mean_p, ce_mean_ps, ce_mean_s, ce_mean_sn, ce_mean_n] = interpolateElectrolyteConcentration(ce,param)
%	interpolateElectrolyteConcentration interpolates the value of electrolyte concentration at the edges of control volumes using harmonic mean.

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

%% Electrolyte concentration interpolation

% Interpolation within the positive electrode
beta_ce_p = 0.5;
ce_mean_p = ce(1:param.Np-1).*ce(2:param.Np)./ (beta_ce_p*ce(2:param.Np) + (1-beta_ce_p)*ce(1:param.Np-1));

% Interpolation on the interface between separator and positive electrode
beta_ce_ps = param.deltax_p*param.len_p/2 / (param.deltax_p*param.len_p/2 + param.deltax_s*param.len_s/2);
ce_mean_ps = ce(param.Np)*ce(param.Np+1)/(beta_ce_ps*ce(param.Np+1) + (1-beta_ce_ps)*ce(param.Np));

% Interpolation within the separator
beta_ce_s = 0.5;
ce_mean_s = ce(param.Np+1:param.Np+param.Ns-1).*ce(param.Np+2:param.Np+param.Ns)./ (beta_ce_s*ce(param.Np+2:param.Np+param.Ns) + (1-beta_ce_s)*ce(param.Np+1:param.Np+param.Ns-1));

% Interpolation on the interface between separator and negative electrode
beta_ce_sn = param.deltax_s*param.len_s/2 / (param.deltax_n*param.len_n/2 + param.deltax_s*param.len_s/2);
ce_mean_sn = ce(param.Np+param.Ns)*ce(param.Np+param.Ns+1)/(beta_ce_sn*ce(param.Np+param.Ns+1) + (1-beta_ce_sn)*ce(param.Np+param.Ns));

% Interpolation within the negative electrode
beta_ce_n = 0.5;
ce_mean_n = ce(param.Np+param.Ns+1:end-1).*ce(param.Np+param.Ns+2:end)./ (beta_ce_n*ce(param.Np+param.Ns+2:end) + (1-beta_ce_n)*ce(param.Np+param.Ns+1:end-1));

end