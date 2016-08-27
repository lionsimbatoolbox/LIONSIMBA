%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IONICFLUX computes the value of the ionic flux in the electrodes. [mol /
% (m^2*s)]
function [jflux,U_p,U_n,dudt_p,dudt_n,J_s] = ionicFlux(ce,cs_star,Phis,Phie,T,solverFlux,film,param)

%% Positive electrode
% Compute the OCV for the positive and negative electrodes.
% [U_p,dudt_p,U_n,dudt_n] = param.OpenCircuitPotentialFunction(cs_star,T,param);
[U_p,dudt_p,U_n,dudt_n] = openCircuitPotential(cs_star,T,param);
% Compute the reaction rates.
[k_pT, k_nT] = param.ReactionRatesFunction(T,param);

% Positive electrode ion flux
deltap = ((0.5*param.F)./(param.R*T(param.Nal+1:param.Nal+param.Np))).*(Phis(1:param.Np)-Phie(1:param.Np)-U_p);
ip = 2*k_pT.*sqrt(ce(1:param.Np)).*sqrt(cs_star(1:param.Np)).*sqrt(param.cs_max(1)-cs_star(1:param.Np));
jnp_calc = ip.* sinh(deltap);

%% Negative electrode

% If ageing is enabled, take into account the SEI resistance
if(param.EnableAgeing==1)
    eta_n   = (Phis(param.Np+1:end)-Phie(param.Np+param.Ns+1:end)-U_n -param.F*solverFlux(param.Np+1:end).*(param.R_SEI+film./(param.k_n_aging)));
else
    eta_n   = (Phis(param.Np+1:end)-Phie(param.Np+param.Ns+1:end)-U_n);
end

deltan      = ((0.5*param.F)./(param.R*T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn))).*eta_n;
in          = 2*k_nT.*sqrt(ce(param.Np+param.Ns+1:end)).*sqrt(cs_star(param.Np+1:end)).*sqrt(param.cs_max(3)-cs_star(param.Np+1:end));
jnn_calc    = in.* sinh(deltan);

J_s = zeros(param.Nn,1);
if(param.EnableAgeing == 1 && param.I > 0)
    eta_s = Phis(param.Np+1:end) - Phie(param.Np+param.Ns+1:end) - param.Uref_s - param.F*solverFlux(param.Np+1:end).*(param.R_SEI+film./(param.k_n_aging));
    % Tafel equation for the side reaction flux. 
    alpha   = 0.5*param.F./(param.R*T(param.Nal+param.Np+param.Ns+1:end-param.Nco));
    J_s = -param.i_0_jside.*(param.I/param.I1C)^param.w.*(exp(-alpha.*eta_s))./param.F;
end
%% Return value
jflux = [jnp_calc;jnn_calc];