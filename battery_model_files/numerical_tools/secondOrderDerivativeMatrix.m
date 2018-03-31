function [derivativeMatrix,r12dxs,dx]=secondOrderDerivativeMatrix(xl,xu,n)
% secondOrderDerivativeMatrix precomputes the matrices used to implement numerical differentiation using 6 points.

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

% Build the overall matrix. The multiplication by 6 is made to avoid
% numerical issues due to the rational values that are present in the
% original matrix definition. A division by 6 is made back in 
% FDM9orderElectrodeDiffusion.m to obtain a proper second order derivative
derivativeMatrix = 6*[first_block;mid_block;last_block];
end