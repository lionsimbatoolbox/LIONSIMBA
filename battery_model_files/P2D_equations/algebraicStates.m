%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGEBRAICSTATES evaluates the residuals of the algebraic equations of the 
% space-discretized version of the P2D model

function [res,Up,Un,dudt_p,dudt_n,Keff,J_S] = algebraicStates(x,ce,cs_barrato,Q,T,film,param,t)

% Ionic flux
jflux   = x(param.jflux_indices-param.ndiff);
% Solid Potential
Phis    = x(param.Phis_indices-param.ndiff);
% Electrolyte potential
Phie    = x(param.Phie_indices-param.ndiff);
% Side reaction flux
js      = x(param.js_indices-param.ndiff);
% Applied current density
I_app   = x(param.Iapp_indices-param.ndiff);
% Check for potentiostatic charge
if(param.AppliedCurrent==3)
    param.I = I_app;
end

% Residuals on the solid potential
res_Phis = solidPhasePotential(jflux + [zeros(param.Np,1);js],...
    param,...
    Phis);

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
    param);

%% Builed the residuals array
% ionic flux residuals
jflux_res = jflux-jflux_calc;
% side reaction residuals
js_res    = js-J_S;

% Check if the potentiostatic mode has been selected.
if(param.AppliedCurrent == 3)
    I_res = Phis(1)-Phis(end)-param.V_reference;
else
    I_res = I_app;
end

% Return the residuals
res = [res_Phie;res_Phis;jflux_res;js_res;I_res];

end
