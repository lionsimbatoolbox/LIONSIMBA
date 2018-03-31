function [D,x] = cheb(N)
%   [D,x] = cheb(N) computes the differentiation matrices for Chebyshev
%   spectral method. Script taken from 
%   https://people.maths.ox.ac.uk/trefethen/spectral.html

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


% computesthe chebyshev differentiation matrix for N points. Also returns node locations
    if N==0, D=0; x=1; return, end
    x = cos(pi*(0:N)/N)';
    c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
    X = repmat(x,1,N+1);
    dX = X-X';
    D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
    D  = D - diag(sum(D'));                 % diagonal entries
end