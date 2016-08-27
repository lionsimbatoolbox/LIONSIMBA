%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLATEELECTROLYTECONCENTRATION evaluates the interpolation of the
% electrolyte concentration at the edges of the control volumes. The
% returned values are obtained using the harmonic mean.

function [ce_mean_p, ce_mean_ps, ce_mean_s, ce_mean_sn, ce_mean_n] = interpolateElectrolyteConcentration(ce,param)
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