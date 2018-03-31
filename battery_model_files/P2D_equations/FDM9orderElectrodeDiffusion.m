function [rhsCs_p, rhsCs_n, ddCs_p, ddCs_n] = FDM9orderElectrodeDiffusion(T, cs_barrato, jflux, dCs, param)
%   FDM9orderElectrodeDiffusion computes the equations for the solid phase
%   diffusion using a Finite Difference scheme with 9 points.ù

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

% This check is performed to ensure a correct sizing of the indices used in the computation of
% the analytical version of the Jacobian matrix.
if(isfield(param,'symbolic_param_num') && (isa(param.Rp_p,'casadi.MX') || isa(param.Rp_p,'casadi.SX') || isa(param.Rp_n,'casadi.MX') || isa(param.Rp_n,'casadi.SX')))
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
ddCs_p 	= [];
ddCs_n 	= [];
start_cs = 1;
% For every single CV in the cathode, let assume the presence of a
% solid particle
for i=1:param.Np
    cs_solid    = cs_barrato((i-1)*param.Nr_p+start_cs:i*param.Nr_p);
    cst_solid   = dCs((i-1)*param.Nr_p+start_cs:i*param.Nr_p);
    
    % Evaluate first order derivatives
    cs_barratox_p                   = param.FO_D_p*cs_solid*param.FO_D_c_p;
    
    % Impose the BCs (the radial direction is normalized between 0 and 1 for improving the numerical robustness)
    % r = 1
    cs_barratox_p(param.Nr_p)       = -jflux(i)/Dps_eff(i)*param.Rp_p;
    
    % r = 0 using l'Hopital rule
    cs_barratox_p(1)                = 0;
    
    % Second order derivatives. The division by 6 comes from the
    % multiplication by 6 of the original differentiation matrix. This
    % reduces the presence of numerical errors.
    cs_barratoxx_p                  = param.SO_D_p*cs_solid*param.SO_D_c_p/6;
    
    % In both r = 0 and r = 1 Neumann BCs are required. For r=0, this
    % process is not carried out since it is required that cs'(r=0) = 0.
    % For r=1 the following modification is performed (according to the
    % particular numerical scheme here adopted for the approximation of the
    % second order derivative)
    cs_barratoxx_p(end)             = cs_barratoxx_p(end) + 50*param.SO_D_dx_p*cs_barratox_p(param.Nr_p)*param.SO_D_c_p;
    
    % Create rhs arrays and residuals
    rhsCs_p_temp    = (Dps_eff(i)*(3*cs_barratoxx_p(1)));
    ddCs_p          = [ddCs_p;param.Rp_p^2*cst_solid(1)-rhsCs_p_temp];
    rhsCs_p         = [rhsCs_p;rhsCs_p_temp];
    rhsCs_p_temp    = (Dps_eff(i)*((cs_barratoxx_p(2:end)+2./(linspace(1/(param.Nr_p-1),1,param.Nr_p-1)').*cs_barratox_p(2:end))));
    ddCs_p          = [ddCs_p;param.Rp_p^2*cst_solid(2:end)-rhsCs_p_temp];
    rhsCs_p         = [rhsCs_p;rhsCs_p_temp];
end

start_cs = param.Np*param.Nr_p+1;

for i=1:param.Nn
    cs_solid    = cs_barrato((i-1)*param.Nr_n+start_cs:i*param.Nr_n+start_cs-1);
    cst_solid   = dCs((i-1)*param.Nr_n+start_cs:i*param.Nr_n+start_cs-1);
    
    % Evaluate first order derivatives
    cs_barratox_n                   = param.FO_D_n*cs_solid*param.FO_D_c_n;
    
    % Impose the BCs (the radial direction is normalized between 0 and 1 for improving the numerical robustness)
    % r = 1
    cs_barratox_n(param.Nr_n)       = -jflux(param.Np+i)/Dns_eff(i)*param.Rp_n;
    
    % r = 0 using l'Hopital rule
    cs_barratox_n(1)                = 0;
    
    % Second order derivatives. The division by 6 comes from the
    % multiplication by 6 of the original differentiation matrix. This
    % reduces the presence of numerical errors.
    cs_barratoxx_n                  = param.SO_D_n*cs_solid*param.SO_D_c_n/6;
    
    % In both r = 0 and r = 1 Neumann BCs are required. For r=0, this
    % process is not carried out since it is required that cs'(r=0) = 0.
    % For r=1 the following modification is performed (according to the
    % particular numerical scheme here adopted for the approximation of the
    % second order derivative)
    cs_barratoxx_n(end)             = cs_barratoxx_n(end) + 50*param.SO_D_dx_n*cs_barratox_n(param.Nr_n)*param.SO_D_c_n;
    
    % Create rhs arrays and residuals
    rhsCs_n_temp    = (Dns_eff(i)*(3*cs_barratoxx_n(1)));
    ddCs_n          = [ddCs_n;param.Rp_n^2*cst_solid(1)-rhsCs_n_temp];
    rhsCs_n         = [rhsCs_n;rhsCs_n_temp];
    rhsCs_n_temp    = (Dns_eff(i)*((cs_barratoxx_n(2:end)+2./(linspace(1/(param.Nr_n-1),1,(param.Nr_n-1))').*cs_barratox_n(2:end))));
    ddCs_n          = [ddCs_n;param.Rp_n^2*cst_solid(2:end)-rhsCs_n_temp];
    rhsCs_n         = [rhsCs_n;rhsCs_n_temp];
end
end