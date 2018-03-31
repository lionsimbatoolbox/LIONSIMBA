function param = solidPhaseDiffusionDifferentiationMatrices(param)
%	solidPhaseDiffusionDifferentiationMatrices pre-computes the numerical scheme for solving the solid phase diffusion equations.
%	These data will be used when Fick's law of diffusion is considered.

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

param.Rad_position_p  = linspace(0,param.Rp_p,param.Nr_p);
if(size(param.Rad_position_p,1)<size(param.Rad_position_p,2))
    param.Rad_position_p = param.Rad_position_p';
end

param.Rad_position_n  = linspace(0,param.Rp_n,param.Nr_n);
if(size(param.Rad_position_n,1)<size(param.Rad_position_n,2))
    param.Rad_position_n = param.Rad_position_n';
end

% Pre-compute the matrices used to compute the first-order numerical
% derivative. For numerical reasons, the radius domain of the particles [0,R] is
% normalized in [0,1]. The numerical formulation is implemented
% accordingly.
[param.FO_D_p,param.FO_D_c_p] = firstOrderDerivativeMatrix(0,1,param.Nr_p);
[param.FO_D_n,param.FO_D_c_n] = firstOrderDerivativeMatrix(0,1,param.Nr_n);

% Pre-compute the matrices used to compute the second-order numerical
% derivative. For numerical reasons, the radius domain of the particles [0,R] is
% normalized in [0,1]. The numerical formulation is implemented
% accordingly.
[param.SO_D_p,param.SO_D_c_p,param.SO_D_dx_p] = secondOrderDerivativeMatrix(0,1,param.Nr_p);
[param.SO_D_n,param.SO_D_c_n,param.SO_D_dx_n] = secondOrderDerivativeMatrix(0,1,param.Nr_n);

end