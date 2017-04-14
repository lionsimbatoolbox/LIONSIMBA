%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECONDORDERDERIVATIVEMATRIX precomputes the matrices used to implement numerical differentiation using 6 points.

function [derivativeMatrix,r12dxs,dx]=secondOrderDerivativeMatrix(xl,xu,n)
%  Grid spacing
dx 		= (xu-xl)/(n-1);
%
r12dxs 	= 1./(12.0*dx^2);

% Define blocks used for the building of the numerical differentiation matrix.
mid_block        = zeros(n-4,n);

%% Numerical scheme
first_row       = [-415/6   +96     -36     +32/3       -3/2    0];
second_row      = [+10      -15     -4      +14         -6      +1];
i_th_row        = [-1       +16     -30     +16         -1];
semi_last_row   = [+1       -6      +14     -4          -15     +10];
last_row        = [0    -3/2    +32/3   -36     +96     -415/6];

%% Blocks building
first_block = [first_row;second_row];
first_block = [first_block zeros(2,n-6)];


last_block 	= [semi_last_row;last_row];
last_block 	= [zeros(2,n-6) last_block];

row_index 	= 1;
for i=3:n-2
    mid_block(row_index,row_index:row_index+4) 	= i_th_row;
    row_index 									= row_index+1;
end

% Build the overall matrix
derivativeMatrix = [first_block;mid_block;last_block];
end