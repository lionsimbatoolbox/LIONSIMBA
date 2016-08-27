%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEATGENERATIONRATES evaluates the heat source terms used in the thermal model.

function [Qrev, Qrxn, Qohm]=heatGenerationRates(Phis,Phie,jflux,T,cs_star,ce,param)

% Retreive effective electrolyte conductivity coefficients.
Keff_p = param.ElectrolyteConductivityFunction(ce(1:param.Np),T(param.Nal+1:param.Nal+param.Np),param,'p');
Keff_s = param.ElectrolyteConductivityFunction(ce(param.Np+1:param.Np+param.Ns),T(param.Nal+param.Np+1:param.Nal+param.Np+param.Ns),param,'s');
Keff_n = param.ElectrolyteConductivityFunction(ce(param.Np+param.Ns+1:end),T(param.Nal+param.Np+param.Ns+1:end-param.Nco),param,'n');


% Evaluate the derivatives used in Qohm calculations
[dPhis, dPhie, dCe] = ThermalDerivatives(Phis',Phie',ce',param);
dPhis = dPhis';
dPhie = dPhie';
dCe = dCe';
[Up,dudt_p,Un,dudt_n] = param.OpenCircuitPotentialFunction(cs_star,T,param);

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
