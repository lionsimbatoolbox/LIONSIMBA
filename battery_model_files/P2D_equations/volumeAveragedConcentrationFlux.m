function [resdQ, rhsQ] = volumeAveragedConcentrationFlux(dQ,Q,jflux,T,param)
%   volumeAveragedConcentrationFlux is used to implement the three parameters reduced model for solid phase diffusion.

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

% Diffusion coefficients for the solid phase
[Dps_eff, Dns_eff] = param.SolidDiffusionCoefficientsFunction(T,param);

rhsQ_p      = -30*Dps_eff./param.Rp_p^2.*Q(1:param.Np)    -45./(2*param.Rp_p^2).*jflux(1:param.Np);

resdQ_p     = dQ(1:param.Np)        -rhsQ_p;

rhsQ_n      = -30*Dns_eff./param.Rp_n^2.*Q(param.Np+1:end)-45./(2*param.Rp_n^2).*jflux(param.Np+1:end);

resdQ_n     = dQ(param.Np+1:end) - rhsQ_n   ;

rhsQ        = [rhsQ_p;rhsQ_n];
resdQ       = [resdQ_p;resdQ_n];

% This model has been taken from the paper, "Efficient Macro-Micro Scale Coupled
% Modeling of Batteries" -  Subramanian,Diwakar,Tapriyal - 2005 JES

end