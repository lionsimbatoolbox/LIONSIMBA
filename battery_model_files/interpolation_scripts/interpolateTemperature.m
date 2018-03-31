function [T_mean_p, T_mean_ps, T_mean_s, T_mean_sn, T_mean_n] = interpolateTemperature(T,param)
%	interpolateTemperature evaluates the interpolation of the temperature at the edges of control volumes using harmonic mean.

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

% Interpolation within the positive electrode
beta_T_p = 0.5;
T_mean_p = T(param.Nal+1:param.Nal+param.Np-1).*T(param.Nal+2:param.Nal+param.Np)./ (beta_T_p*T(param.Nal+2:param.Nal+param.Np) + (1-beta_T_p)*T(param.Nal+1:param.Nal+param.Np-1));

% Interpolation on the interface between separator and positive electrode
beta_T_ps = param.deltax_p*param.len_p/2 / (param.deltax_p*param.len_p/2 + param.deltax_s*param.len_s/2);
T_mean_ps = T(param.Nal+param.Np)*T(param.Nal+param.Np+1)/(beta_T_ps*T(param.Nal+param.Np+1) + (1-beta_T_ps)*T(param.Nal+param.Np));

% Interpolation within the separator
beta_T_s = 0.5;
T_mean_s = T(param.Nal+param.Np+1:param.Nal+param.Np+param.Ns-1).*T(param.Nal+param.Np+2:param.Nal+param.Np+param.Ns)./ (beta_T_s*T(param.Nal+param.Np+2:param.Nal+param.Np+param.Ns) + (1-beta_T_s)*T(param.Nal+param.Np+1:param.Nal+param.Np+param.Ns-1));

% Interpolation on the interface between separator and negative electrode
beta_T_sn = param.deltax_s*param.len_s/2 / (param.deltax_n*param.len_n/2 + param.deltax_s*param.len_s/2);
T_mean_sn = T(param.Nal+param.Np+param.Ns)*T(param.Nal+param.Np+param.Ns+1)/(beta_T_sn*T(param.Nal+param.Np+param.Ns+1) + (1-beta_T_sn)*T(param.Nal+param.Np+param.Ns));

% Interpolation within the negative electrode
beta_T_n = 0.5;
T_mean_n = T(param.Nal+param.Np+param.Ns+1:end-(param.Ncu)-1).*T(param.Nal+param.Np+param.Ns+2:end-param.Ncu)./ (beta_T_n*T(param.Nal+param.Np+param.Ns+2:end-param.Ncu) + (1-beta_T_n)*T(param.Nal+param.Np+param.Ns+1:end-(param.Ncu+1)));

end