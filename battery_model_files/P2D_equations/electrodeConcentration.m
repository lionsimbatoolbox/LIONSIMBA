%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELECTRODECONCENTRATION evaluates the solid phase concentration of Li-ions.

function [ddCs, rhsCs] = electrodeConcentration(dCs,cs_barrato,T,jflux,param)
% ODE descsribing the concentration of lithium ions within the
% electrodes.

if(param.SolidPhaseDiffusion==1 || param.SolidPhaseDiffusion==2)
    % Cathode
    rhsCs_p  =((-3/param.Rp_p)*jflux(1:param.Np));
    ddCs_p   = dCs(1:param.Np) - rhsCs_p;
    
    % Anode
    rhsCs_n  = ((-3/param.Rp_n)*jflux(param.Np+1:end));
    ddCs_n   = dCs(param.Np+1:end) - rhsCs_n;
else
	if(isfield(param,'symbolic_param_num') && (isa(param.Rp_p,'casadi.MX') || isa(param.Rp_p,'casadi.SX') || isa(param.Rp_n,'casadi.MX') || isa(param.Rp_n,'casadi.SX')))
        % If the simulation is used for sensitivity purposes, the radius
        % positions need to be evaluated at each time.
        param.Rad_position_p  = linspace(0,param.Rp_p,param.Nr_p);
        param.Rad_position_n  = linspace(0,param.Rp_n,param.Nr_n);
    end
    % If the regular diffusion equation is selected, then use FDM to
    % evaluate the complete solution.
    % First, retreive the diffusion coefficients
    [Dps_eff, Dns_eff] = param.SolidDiffusionCoefficientsFunction(T,param);
    % Initialize the variables
    rhsCs_p = [];
    rhsCs_n = [];
    ddCs_p = [];
    ddCs_n = [];
    start_cs = 1;
    % For every single CV in the cathode, let assume the presence of a
    % solid particle
    for i=1:param.Np
        cs_solid    = cs_barrato((i-1)*param.Nr_p+start_cs:i*param.Nr_p);
        cst_solid   = dCs((i-1)*param.Nr_p+start_cs:i*param.Nr_p);
        % Evaluate first order derivatives
        %         cs_barratox_p                   = deriveFirst(0,param.Rp_p,param.Nr_p,cs_solid)';
        cs_barratox_p                   = param.FO_D_p*cs_solid*param.FO_D_c_p;
        % Impose the BCs
        % r=Rp_p
        cs_barratox_p(param.Nr_p)       = -jflux(i)/Dps_eff(i);
        % r=0 using l'Hopitàl rule
        cs_barratox_p(1)                = 0;
        % Second order derivatives
        cs_barratoxx_p                  = param.SO_D_p*cs_solid*param.SO_D_c_p;
        % According to the numerical scheme for the second order derivative, add the term to enforce Neumann BCs. Note that another
        % Neumann BC should be imposed at cs_barratoxx_p but given that at that point the derivative is zero, its contribution here
        % is neglected.
        cs_barratoxx_p(end)             = cs_barratoxx_p(end) + 50*param.SO_D_dx_p*cs_barratox_p(param.Nr_p)*param.SO_D_c_p;
        rhsCs_p_temp    = (Dps_eff(i)*(3*cs_barratoxx_p(1)));
        ddCs_p          = [ddCs_p;cst_solid(1)-rhsCs_p_temp];
        rhsCs_p         = [rhsCs_p;rhsCs_p_temp];
        rhsCs_p_temp    = (Dps_eff(i)*((cs_barratoxx_p(2:end)+2./(param.Rad_position_p(2:end)).*cs_barratox_p(2:end))));
        ddCs_p          = [ddCs_p;cst_solid(2:end)-rhsCs_p_temp];
        rhsCs_p         = [rhsCs_p;rhsCs_p_temp];
    end
    
    start_cs = param.Np*param.Nr_p+1;
    
    for i=1:param.Nn
        cs_solid    = cs_barrato((i-1)*param.Nr_n+start_cs:i*param.Nr_n+start_cs-1);
        cst_solid   = dCs((i-1)*param.Nr_n+start_cs:i*param.Nr_n+start_cs-1);
        % Evaluate first order derivatives
        %         cs_barratox_n                   = deriveFirst(0,param.Rp_n,param.Nr_n,cs_solid)';
        cs_barratox_n                   = param.FO_D_n*cs_solid*param.FO_D_c_n;
        % Impose the BCs
        % r=Rp_p
        cs_barratox_n(param.Nr_n)       = -jflux(param.Np+i)/Dns_eff(i);
        % r=0 using l'Hopitàl rule
        cs_barratox_n(1)                = 0;
        % According to the numerical scheme for the second order derivative, add the term to enforce Neumann BCs. Note that another
        % Neumann BC should be imposed at cs_barratoxx_p but given that at that point the derivative is zero, its contribution here
        % is neglected.
        cs_barratoxx_n                  = param.SO_D_n*cs_solid*param.SO_D_c_n;
        cs_barratoxx_n(end)             = cs_barratoxx_n(end) + 50*param.SO_D_dx_n*cs_barratox_n(param.Nr_n)*param.SO_D_c_n;
        rhsCs_n_temp    = (Dns_eff(i)*(3*cs_barratoxx_n(1)));
        ddCs_n          = [ddCs_n;cst_solid(1)-rhsCs_n_temp];
        rhsCs_n         = [rhsCs_n;rhsCs_n_temp];
        rhsCs_n_temp    = (Dns_eff(i)*((cs_barratoxx_n(2:end)+2./(param.Rad_position_n(2:end)).*cs_barratox_n(2:end))));
        ddCs_n          = [ddCs_n;cst_solid(2:end)-rhsCs_n_temp];
        rhsCs_n         = [rhsCs_n;rhsCs_n_temp];
    end
end
rhsCs   = [rhsCs_p;rhsCs_n];
ddCs    = [ddCs_p;ddCs_n];
end