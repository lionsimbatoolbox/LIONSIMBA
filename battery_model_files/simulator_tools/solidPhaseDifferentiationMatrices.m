%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function param = solidPhaseDifferentiationMatrices(param)
% Phis matrices preallocation. Given that in the proposed code the
% solid phase conductivity is considered to be constant across the cell
% length, the discretization matrices can be assembled only once and
% used later in the code.

% A matrix for the positive electrode
c = ones(param.Np-1,1);
d = -2*ones(param.Np,1);

A_p 				= gallery('tridiag',c,d,c);
A_p(1,1) 			= -1;
A_p(end,end-1:end) 	= [1 -1];

% A matrix for the negative electrode
c = ones(param.Nn-1,1);
d = -2*ones(param.Nn,1);

A_n 			= gallery('tridiag',c,d,c);
A_n(1,1) 		= -1;
A_n(end,end) 	= -1;

% Store the matrices in the param structure for future usage for each
% cell in the pack.
param.A_p = A_p;
param.A_n = A_n;
end