%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function param = solidPhaseDiffusionDifferentiationMatrices(param)

% Precompute the discretization points for the solid particles.
% These data will be used when Fick's law of diffusion is considered.
param.Rad_position_p  = linspace(0,param.Rp_p,param.Nr_p);
if(size(param.Rad_position_p,1)<size(param.Rad_position_p,2))
    param.Rad_position_p = param.Rad_position_p';
end

param.Rad_position_n  = linspace(0,param.Rp_n,param.Nr_n);
if(size(param.Rad_position_n,1)<size(param.Rad_position_n,2))
    param.Rad_position_n = param.Rad_position_n';
end

% Precompute the matrices used for the numerical differentiation. These matrices
% will be used when Fick's law is considered.
[param.FO_D_p,param.FO_D_c_p] = firstOrderDerivativeMatrix(0,param.Rp_p,param.Nr_p);
[param.FO_D_n,param.FO_D_c_n] = firstOrderDerivativeMatrix(0,param.Rp_n,param.Nr_n);

% Precompute the matrices used for the numerical differentiation. These matrices
% will be used when Fick's law is considered.
[param.SO_D_p,param.SO_D_c_p,param.SO_D_dx_p] = secondOrderDerivativeMatrix(0,param.Rp_p,param.Nr_p);
[param.SO_D_n,param.SO_D_c_n,param.SO_D_dx_n] = secondOrderDerivativeMatrix(0,param.Rp_n,param.Nr_n);

end