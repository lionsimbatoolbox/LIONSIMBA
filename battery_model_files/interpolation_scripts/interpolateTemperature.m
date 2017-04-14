%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLATETEMPERATURE evaluates the interpolation of the temperature at
% the edges of the control volumes.

function [T_mean_p, T_mean_ps, T_mean_s, T_mean_sn, T_mean_n] = interpolateTemperature(T,param)

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
T_mean_n = T(param.Nal+param.Np+param.Ns+1:end-(param.Nco)-1).*T(param.Nal+param.Np+param.Ns+2:end-param.Nco)./ (beta_T_n*T(param.Nal+param.Np+param.Ns+2:end-param.Nco) + (1-beta_T_n)*T(param.Nal+param.Np+param.Ns+1:end-(param.Nco+1)));

end