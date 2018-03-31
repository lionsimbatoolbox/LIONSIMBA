function [jflux,U_p,U_n,dudt_p,dudt_n,J_s] = ionicFlux(ce,cs_star,Phis,Phie,T,solverFlux,film,param,sign_input_density,I_density)
% ionicFlux Computes the molar flux density of Li-ions at the electrode-electrolyte interface [mol/(m^2*s)].

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

%% Positive electrode
% Compute the OCV for the positive and negative electrodes.
[U_p,dudt_p,U_n,dudt_n] = param.OpenCircuitPotentialFunction(cs_star,T,param,sign_input_density);

% Compute the reaction rates.
[k_pT, k_nT] = param.ReactionRatesFunction(T,param);

% Positive electrode ion flux
deltap = ((0.5*param.F)./(param.R*T(param.Nal+1:param.Nal+param.Np))).*(Phis(1:param.Np)-Phie(1:param.Np)-U_p);
ip = 2*k_pT.*sqrt(ce(1:param.Np)).*sqrt(cs_star(1:param.Np)).*sqrt(param.cs_maxp-cs_star(1:param.Np));
jnp_calc = ip.* sinh(deltap);

%% Negative electrode

% If ageing is enabled, take into account the SEI resistance
if(param.EnableAgeing==1)
    eta_n   = (Phis(param.Np+1:end)-Phie(param.Np+param.Ns+1:end)-U_n -param.F*solverFlux(param.Np+1:end).*(param.R_SEI+film./(param.k_n_aging)));
else
    eta_n   = (Phis(param.Np+1:end)-Phie(param.Np+param.Ns+1:end)-U_n);
end

deltan      = ((0.5*param.F)./(param.R*T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn))).*eta_n;
in          = 2*k_nT.*sqrt(ce(param.Np+param.Ns+1:end)).*sqrt(cs_star(param.Np+1:end)).*sqrt(param.cs_maxn-cs_star(param.Np+1:end));
jnn_calc    = in.* sinh(deltan);

J_s = zeros(param.Nn,1);

% Switch cases when the applied current density is a symbolic variable
if(isa(I_density,'casadi.MX')||isa(I_density,'casadi.SX') && param.EnableAgeing == 1)
    eta_s = Phis(param.Np+1:end) - Phie(param.Np+param.Ns+1:end) - param.Uref_s - param.F*solverFlux(param.Np+1:end).*(param.R_SEI+film./(param.k_n_aging));
    % Tafel equation for the side reaction flux.
    alpha   = 0.5*param.F./(param.R*T(param.Nal+param.Np+param.Ns+1:end-param.Ncu));
    % By means of the if_else statement of CasADi, it is possible to represent dynamics that switch according to the value of the symbolic quantity I_density
    J_s = if_else(I_density>=0,-param.i_0_jside.*(I_density/param.I1C)^param.w.*(exp(-alpha.*eta_s))./param.F,zeros(param.Nn,1));
elseif(param.EnableAgeing == 1 && sign_input_density > 0)
    eta_s = Phis(param.Np+1:end) - Phie(param.Np+param.Ns+1:end) - param.Uref_s - param.F*solverFlux(param.Np+1:end).*(param.R_SEI+film./(param.k_n_aging));
    % Tafel equation for the side reaction flux.
    alpha   = 0.5*param.F./(param.R*T(param.Nal+param.Np+param.Ns+1:end-param.Ncu));
    J_s = -param.i_0_jside.*(I_density/param.I1C)^param.w.*(exp(-alpha.*eta_s))./param.F;
end
%% Return value
jflux = [jnp_calc;jnn_calc];
