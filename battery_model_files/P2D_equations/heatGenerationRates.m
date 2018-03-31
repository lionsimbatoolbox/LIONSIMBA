function [Qrev, Qrxn, Qohm]=heatGenerationRates(Phis,Phie,jflux,T,cs_star,ce,param)
% heatGenerationRates evaluates the heat source terms used in the thermal model.

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

% Retrieve effective electrolyte conductivity coefficients.
Keff_p = param.ElectrolyteConductivityFunction(ce(1:param.Np),T(param.Nal+1:param.Nal+param.Np),param,'p');
Keff_s = param.ElectrolyteConductivityFunction(ce(param.Np+1:param.Np+param.Ns),T(param.Nal+param.Np+1:param.Nal+param.Np+param.Ns),param,'s');
Keff_n = param.ElectrolyteConductivityFunction(ce(param.Np+param.Ns+1:end),T(param.Nal+param.Np+param.Ns+1:end-param.Ncu),param,'n');

% Evaluate the derivatives used in Qohm calculations
[dPhis, dPhie, dCe] = ThermalDerivatives(Phis',Phie',ce',param);	% y(end) obtains current from a vector containing algebraic states
dPhis = dPhis';
dPhie = dPhie';
dCe = dCe';
[Up,dudt_p,Un,dudt_n] = param.OpenCircuitPotentialFunction(cs_star,T,param,sign_input_density);

%% Reversible heat generation rate

% Positive electrode
Qrev_p = param.F*param.a_i(1)*jflux(1:param.Np).*T(param.Nal+1:param.Nal+param.Np).*dudt_p;

% Negative Electrode
Qrev_n = param.F*param.a_i(3)*jflux(param.Np+1:end).*T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn).*dudt_n;

%% Reaction heat generation rate

% Positive overpotential
eta_p = (Phis(1:param.Np)-Phie(1:param.Np)-Up);
% Positive reaction heat generation rate
Qrxn_p = param.F*param.a_i(1)*jflux(1:param.Np).*eta_p;

% Negative overpotential
eta_n = (Phis(param.Np+1:end)-Phie(param.Np+param.Ns+1:end)-Un);
% Negative reaction heat generation rate
Qrxn_n = param.F*param.a_i(3)*jflux(param.Np+1:end).*eta_n;

%% Ohmic heat generation rate

% Positive electrode ohmic generation rate
Qohm_p = param.sig_eff(1) * (dPhis(1:param.Np)).^2 + Keff_p.*(dPhie(1:param.Np)).^2 + 2*param.R*Keff_p.*T(param.Nal+1:param.Nal+param.Np)*(1-param.tplus)/param.F.*dCe(1:param.Np).*1./ce(1:param.Np).*dPhie(1:param.Np);
% Separator ohmic generation rate
Qohm_s = Keff_s.*(dPhie(param.Np+1:param.Np+param.Ns)).^2 + 2*param.R*Keff_s.*T(param.Nal+param.Np+1:param.Nal+param.Np+param.Ns)*(1-param.tplus)/param.F.*dCe(param.Np+1:param.Np+param.Ns).*1./ce(param.Np+1:param.Np+param.Ns).*dPhie(param.Np+1:param.Np+param.Ns);
% Negative electrode ohmic generation rate
Qohm_n = param.sig_eff(3) * (dPhis(param.Np+1:end)).^2 +Keff_n.*(dPhie(param.Np+param.Ns+1:end)).^2 + 2*param.R*Keff_n.*T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn)*(1-param.tplus)/param.F.*dCe(param.Np+param.Ns+1:end).*1./ce(param.Np+param.Ns+1:end).*dPhie(param.Np+param.Ns+1:end);

Qrev = [Qrev_p zeros(1,param.Ns) Qrev_n];
Qrxn = [Qrxn_p zeros(1,param.Ns) Qrxn_n];
Qohm = [Qohm_p Qohm_s Qohm_n];
