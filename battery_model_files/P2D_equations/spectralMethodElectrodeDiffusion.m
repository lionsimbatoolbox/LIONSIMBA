function [rhsCs_p, rhsCs_n, ddCs_p, ddCs_n] = spectralMethodElectrodeDiffusion(T, cs_barrato, jflux, dCs, param)
%   spectralMethodElectrodeDiffusion computes the solid phase diffusion
%   using a spectral scheme (not fully completed).

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

% First, retreive the diffusion coefficients
[Dps_eff, Dns_eff] = param.SolidDiffusionCoefficientsFunction(T,param);
% Initialize the variables
rhsCs_p = [];
rhsCs_n = [];
ddCs_p = [];
ddCs_n = [];
start_cs = 1;

[DiffMat_p,param.Rad_position_p] = cheb(param.Nr_p-1);
[DiffMat_n,param.Rad_position_n] = cheb(param.Nr_n-1);
param.Rad_position_p_flipped = flipud(param.Rad_position_p); % flipped to be LIONSIMBA compatible
param.Rad_position_n_flipped = flipud(param.Rad_position_n);

% For every single CV in the cathode, let assume the presence of a solid particle
for i=1:param.Np
    % Lionsimba scheme vector: 1 (top) is particle centre, end(bottom) - particle surface
    cs_solid    = cs_barrato((i-1)*param.Nr_p+start_cs:i*param.Nr_p);
    cst_solid   = dCs((i-1)*param.Nr_p+start_cs:i*param.Nr_p);
    
    cs_solid_flipped   = flipud(cs_solid);      % cheb matrices work on [1 to -1] ordering
    cs_barratox_p      = DiffMat_p*cs_solid_flipped;
    cs_barratox_p(1)   = -jflux(i)*param.Rp_p*0.5/Dps_eff(i); % modified BC value due to cheb scheme
    cs_barratox_p(end) = 0;
    
    % Below line: we compute the RHS of the Fick's law PDE as it appears in LIONSIMBA
    % paper, but in scaled co-ordinate system. Note that this is missing a leading
    % scaling factor (1/(r+1)^2),which will be included in later lines of code in
    % this file At the centre, the equation becomes degenerate and hence treated
    % separately. The scaling factor is included at a later point below.
    
    rhs_numerator_p = DiffMat_p*(4*Dps_eff(i)*((param.Rad_position_p + 1).^2).*cs_barratox_p/(param.Rp_p^2));
    
    rhs_limit_vector = (4*Dps_eff(i)/param.Rp_p^2)*3*(DiffMat_p*cs_barratox_p); % limit at r_tilde tends to -1 (at centre)
    cs_barratox_p  = flipud(cs_barratox_p); % flip back to be compatible with LIONSIMBA
    rhs_numerator_p = flipud(rhs_numerator_p);
    
    rhsCs_p_temp    = rhs_limit_vector(end); % clever trick to apply the L'hopital's rule at particle centre 'end' is used because of cpu time needed for flipud
    ddCs_p          = [ddCs_p;cst_solid(1)-rhsCs_p_temp];
    rhsCs_p         = [rhsCs_p;rhsCs_p_temp];
    rhsCs_p_temp    = rhs_numerator_p(2:end)./((param.Rad_position_p_flipped(2:end) + 1).^2); % Apply the scaling factor for the Fick's law RHS for the rest of the shells excluding particle centre
    ddCs_p          = [ddCs_p;cst_solid(2:end)-rhsCs_p_temp];
    rhsCs_p         = [rhsCs_p;rhsCs_p_temp];
end

start_cs = param.Np*param.Nr_p+1;

for i=1:param.Nn
    % Lionsimba scheme vector: 1 (top) is particle centre, end(bottom) - particle surface
    cs_solid    = cs_barrato((i-1)*param.Nr_n+start_cs:i*param.Nr_n+start_cs-1);
    cst_solid   = dCs((i-1)*param.Nr_n+start_cs:i*param.Nr_n+start_cs-1);
    
    cs_solid_flipped   = flipud(cs_solid);      % cheb matrices work with [1 to -1] ordering
    cs_barratox_n      = DiffMat_n*cs_solid_flipped;
    cs_barratox_n(1)   = -jflux(param.Np+i)*param.Rp_n*0.5/Dns_eff(i); % modified BC value due to cheb scheme
    cs_barratox_n(end) = 0;
    
    % Below line: we compute the RHS of the Fick's law PDE as it appears in LIONSIMBA
    % paper, but in scaled co-ordinate system. Note that this is missing a leading
    % scaling factor (1/(r+1)^2) which will be accounted for in later lines of code
    % in this file, At the centre, the equation becomes degenerate and hence treated
    % separately. The scaling factor is included at a later point below.
    
    rhs_numerator_n = DiffMat_n*(4*Dns_eff(i)*((param.Rad_position_n + 1).^2).*cs_barratox_n/(param.Rp_n^2));
    
    rhs_limit_vector = (4*Dns_eff(i)/param.Rp_n^2)*3*(DiffMat_n*cs_barratox_n); % limit at r_tilde tends to -1 (at centre)
    cs_barratox_n  = flipud(cs_barratox_n); % flip back to be compatible with LIONSIMBA
    rhs_numerator_n = flipud(rhs_numerator_n);
    
    rhsCs_n_temp    = rhs_limit_vector(end); % clever trick to apply the L'hopital's rule at particle centre 'end' is used because of cpu time needed for flipud
    ddCs_n          = [ddCs_n;cst_solid(1)-rhsCs_n_temp];
    rhsCs_n         = [rhsCs_n;rhsCs_n_temp];
    rhsCs_n_temp    = rhs_numerator_n(2:end)./((param.Rad_position_n_flipped(2:end) + 1).^2); % Apply the scaling factor for the Fick's law RHS for the rest of the shells excluding particle centre
    ddCs_n          = [ddCs_n;cst_solid(2:end)-rhsCs_n_temp];
    rhsCs_n         = [rhsCs_n;rhsCs_n_temp];
end
end

