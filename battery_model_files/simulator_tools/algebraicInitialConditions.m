%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x0_alg, n_alg] =   algebraicInitialConditions(param)
% Initial guess for the algebraic variables
jflux_init          = [-0.43e-5*ones(param.Np,1);0.483e-5*ones(param.Nn,1)];
Phis_init           = [4.2*ones(param.Np,1);0.074*ones(param.Nn,1)];
Phie_init           = zeros(param.Nsum,1);
js_init             = 0.483e-5*ones(param.Nn,1);
I_app               = 1;

% Build the array of algebraic initial conditions
x0_alg              = [
    jflux_init;...
    Phis_init;...
    Phie_init;...
    js_init;...
    I_app
    ];

n_alg = length(x0_alg);
end