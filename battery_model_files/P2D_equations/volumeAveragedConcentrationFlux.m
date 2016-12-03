%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% volumeAveragedConcentrationFlux is used to implement the three parameters reduced model for solid phase diffusion.

function [resdQ, rhsQ] = volumeAveragedConcentrationFlux(dQ,Q,jflux,T,param)
% Diffusion coefficients for the solid phase
[Dps_eff, Dns_eff] = param.SolidDiffusionCoefficientsFunction(T,param);

rhsQ_p      = -30*Dps_eff./param.Rp_p^2.*Q(1:param.Np)    -45./(2*param.Rp_p^2).*jflux(1:param.Np);

resdQ_p     = dQ(1:param.Np)        -rhsQ_p;

rhsQ_n      = -30*Dns_eff./param.Rp_n^2.*Q(param.Np+1:end)-45./(2*param.Rp_n^2).*jflux(param.Np+1:end);

resdQ_n     = dQ(param.Np+1:end) - rhsQ_n   ;

rhsQ        = [rhsQ_p;rhsQ_n];
resdQ       = [resdQ_p;resdQ_n];

% This model has been taken from "Efficient Macro-Micro Scale Coupled
% Modeling of Batteries" -  Subramanian,Diwakar,Tapriyal - 2005 JES

end