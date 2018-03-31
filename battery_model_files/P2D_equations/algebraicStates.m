function [res,Up,Un,dudt_p,dudt_n,Keff,J_S] = algebraicStates(x,ce,cs_barrato,Q,T,film,param,t)
% algebraicStates evaluates the residuals of all the algebraic equations at respective cell centres.

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

% Ionic flux
jflux     = x(param.jflux_indices-param.ndiff);
% Solid Potential
Phis      = x(param.Phis_indices-param.ndiff);
% Electrolyte potential
Phie      = x(param.Phie_indices-param.ndiff);
% Side reaction flux
js        = x(param.js_indices-param.ndiff);

I_density = x(param.curr_dens_indices-param.ndiff);

sign_input_density = evaluate_sign_input_density(param); % evaluates sign of input current/power density as per operating mode

if param.OperatingMode==3 % Check if operating in potentiostatic charge mode
    param.I_density = I_density;
end

% Residuals on the solid potential
res_Phis = solidPhasePotential(jflux + [zeros(param.Np,1);js],...
    param,...
    Phis,I_density);

% Residuals on the electrolyte potential
[res_Phie, Keff] = electrolytePotential(jflux + [zeros(param.Np,1);js],...
    ce,...
    T,...
    param,...
    Phie);

% Surface average concentration
cs_star = surfaceConcentration(cs_barrato,...
    jflux,...
    Q,...
    T,...
    param);

% Ionic flux calculations
[jflux_calc,Up,Un,dudt_p,dudt_n,J_S] = ionicFlux(ce,...
    cs_star,...
    Phis,...
    Phie,...
    T,...
    jflux + [zeros(param.Np,1);js],...
    film,...
    param,sign_input_density,I_density);

%% Build the residuals array
% ionic flux residuals
jflux_res = jflux-jflux_calc;
% side reaction residuals
js_res    = js-J_S;


switch param.edge_values
	% By linear interpolation
	case 2
		phis_pos_cc = 1.5*Phis(1) - 0.5*Phis(2);
		phis_neg_cc = -0.5*Phis(end-1) + 1.5*Phis(end);
	% Consider the centroid value
	otherwise
		phis_pos_cc = Phis(1);
		phis_neg_cc = Phis(end);
end

% The variable 'scalar_res' below is a place-holder for the scalar residual variable that represents either the residual of current density (in modes 1,2,4 and 5) or voltage residual (in mode 3)
if param.OperatingMode==1 || param.OperatingMode==4
    scalar_res = I_density-param.I_density;
    res = [res_Phie;res_Phis;jflux_res;js_res;scalar_res]; % Return the residuals (note that res_Phis includes the boundary conditions for both electrodes as per applied current density)
elseif param.OperatingMode==2 || param.OperatingMode==5
    scalar_res = I_density-(param.P_density/(phis_pos_cc-phis_neg_cc));
    pwr_dens_neg_BC_res = + (param.P_density) ...
                          + (param.sig_eff(1)/(param.deltax_p*param.len_p))*(Phis(2)*phis_pos_cc-Phis(1)*phis_pos_cc) ...
                          + (param.sig_eff(3)/(param.deltax_n*param.len_n))*(Phis(end-1)*phis_neg_cc-Phis(end)*phis_neg_cc) ...
                          - (param.len_p*param.deltax_p*param.a_i(1)*param.F*jflux(1))*phis_pos_cc ...
                          - (param.len_n*param.deltax_n*param.a_i(3)*param.F*jflux(end))*phis_neg_cc;
    res = [res_Phie;res_Phis;pwr_dens_neg_BC_res;jflux_res;js_res;scalar_res]; % Return the residuals (note that the power density residual at the negative electrode has been included in the residual vector)
else
    scalar_res = phis_pos_cc-phis_neg_cc-param.V_reference;
    res = [res_Phie;res_Phis;jflux_res;js_res;scalar_res]; % Need to validate this. Return the Residuals (for the potentiostatic case)
end
end
