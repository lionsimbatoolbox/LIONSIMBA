%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLATEELECTROLYTECONCENTRATIONFLUXES evaluates the interpolation of the
% electrolyte concentration flux at the edges of the control volumes. The
% returned values are obtained using the harmonic mean.

function [ce_flux_p, ce_flux_ps, ce_flux_s, ce_flux_sn, ce_flux_n] = interpolateElectrolyteConcetrationFluxes(ce,param)

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