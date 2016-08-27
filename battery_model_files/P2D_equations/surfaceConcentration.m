%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SURFACECONCENTRATION evaluates the average surface concentration of
% Li-ions within the electrodes.

function cs_star = surfaceConcentration(cs_barrato,jflux,Q,T,param)
% Diffusion coefficients for the solid phase
[Dps_eff, Dns_eff] = param.SolidDiffusionCoefficientsFunction(T,param);

% Check what kind of solid diffusion model has been chosen.
if(param.SolidPhaseDiffusion==1) % Two parameters model
    % Evaluates the average surface concentration in both the electrodes.
    % Cathode
    cs_star_p = cs_barrato(1:param.Np)-(param.Rp_p./(Dps_eff.*5)).*jflux(1:param.Np);
    % Anode
    cs_star_n = cs_barrato(param.Np+1:end)-(param.Rp_n./(Dns_eff.*5)).*jflux(param.Np+1:end);
elseif(param.SolidPhaseDiffusion==2) % Three parameters model
    % Cathode
    cs_star_p = cs_barrato(1:param.Np)+(param.Rp_p./(Dps_eff.*35)).*(-jflux(1:param.Np)+8*Dps_eff.*Q(1:param.Np));
    % Anode
    cs_star_n = cs_barrato(param.Np+1:end)+(param.Rp_n./(Dns_eff.*35)).*(-jflux(param.Np+1:end)+8*Dns_eff.*Q(param.Np+1:end));
elseif(param.SolidPhaseDiffusion==3) % Full model
    p_indices = param.Nr_p:param.Nr_p:param.Nr_p*param.Np;
    n_indices = param.Nr_n:param.Nr_n:param.Nr_n*param.Nn;
    % If the full model has been used, just take the concentration data at
    % r=Rp
    cs_star_p = cs_barrato(p_indices);
    cs_star_n = cs_barrato(n_indices+p_indices(end));
end
% Return the residuals.
cs_star = [cs_star_p;cs_star_n];

end