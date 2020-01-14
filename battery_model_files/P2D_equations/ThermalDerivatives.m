function [dPhis, dPhie, dCe] = ThermalDerivatives(Phis,Phie,ce,param)
% ThermalDerivatives evaluates the set of spatial derivatives used in thermal dynamics.
% This function evaluates the spatial derivatives (gradients) of the relevant variables for
% thermal (heat generation) calculations. Note that the spatial derivatives are evaluated within,
% and not at the edges of control volumes.

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

% For each of the numerical derivatives computed below, the first and last control volumes are evaluated with first
% order accuracy (forward and backward difference schemes respectively),
% while the middle control volume approximations use a second order accuracy (central difference scheme).

%% Solid potential derivatives

% Positive Electrode
dPhisp = [(-3*Phis(1)+4*Phis(2)-Phis(3))/(2*param.deltax_p*param.len_p);... 						% Forward differentiation scheme
    (Phis(3:param.Np)-Phis(1:param.Np-2)) / (2*param.deltax_p*param.len_p);...						% Central differentiation scheme
    (3*Phis(param.Np)-4*Phis(param.Np-1)+Phis(param.Np-2)) / (2*param.deltax_p*param.len_p);...		% Backward differentiation scheme
    ];

% Negative Electrode
dPhisn = [(-3*Phis(param.Np+1)+4*Phis(param.Np+2)-Phis(param.Np+3))/(2*param.deltax_n*param.len_n);... 	% Forward differentiation scheme
    (Phis(param.Np+3:end)-Phis(param.Np+1:end-2)) / (2*param.deltax_n*param.len_n);... 					% Central differentiation scheme
    (3*Phis(end)-4*Phis(end-1)+Phis(end-2)) / (2*param.deltax_n*param.len_n);...						% Backward differentiation scheme
    ];

dPhis = [dPhisp;dPhisn];

%% Electrolyte potential derivatives

% Positive Electrode

dPhiep = [ (-3*Phie(1)+4*Phie(2)-Phie(3))/(2*param.deltax_p*param.len_p);...		% Forward differentiation scheme
    (Phie(3:param.Np)-Phie(1:param.Np-2))/(2*param.deltax_p*param.len_p);...	  	% Central differentiation scheme
    ];

% Attention! The last volume of the positive electrode will involve one volume of the
% separator for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.

% Last CV in the positive electrode: derivative approximation with a central scheme
dPhie_last_p = 2*(Phie(param.Np+1)-Phie(param.Np-1))/(3 * param.deltax_p*param.len_p + param.deltax_s*param.len_s);

% Separator

% Attention! The first volume of the separator will involve one volume of the
% positive section for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.

% First CV in the separator: derivative approximation with a central difference scheme
dPhie_first_s = 2*(Phie(param.Np+2)-Phie(param.Np))/(param.deltax_p*param.len_p + 3* param.deltax_s*param.len_s);

% Central difference scheme
dPhies =  (Phie(param.Np+3:param.Np+param.Ns)-Phie(param.Np+1:param.Np+param.Ns-2))/(2*param.deltax_s*param.len_s);

% Attention! The last volume of the separator will involve one volume of the
% negative section for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.

% Last CV in the separator: derivative approximation with a central scheme
dPhie_last_s = 2*(Phie(param.Np+param.Ns+1)-Phie(param.Np+param.Ns-1))/( param.deltax_n*param.len_n + 3*param.deltax_s*param.len_s);

% Negative electrode

% Attention! The first volume of the negative electrode will involve one volume of the
% separator section for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.

% First CV in the negative electrode: derivative approximation with a central scheme
dPhie_first_n = 2*(Phie(param.Np+param.Ns+2)-Phie(param.Np+param.Ns))/(3 * param.deltax_n*param.len_n + param.deltax_s*param.len_s);

% Central difference scheme
dPhien = [(Phie(param.Np+param.Ns+3:end)-Phie(param.Np+param.Ns+1:end-2))/(2*param.deltax_n*param.len_n);...
    (3*Phie(end)-4*Phie(end-1)+Phie(end-2))/(2*param.deltax_n*param.len_n)
    ];

dPhie = [dPhiep;dPhie_last_p;dPhie_first_s;dPhies;dPhie_last_s;dPhie_first_n;dPhien];

%% Electrolyte concentration derivatives

% Positive Electrode

dCep = [ (-3*ce(1)+4*ce(2)-ce(3))/(2*param.deltax_p*param.len_p);... 		% Forward differentiation scheme
    (ce(3:param.Np)-ce(1:param.Np-2))/(2*param.deltax_p*param.len_p);... 	% Central differentiation scheme
    ];

% Attention! The last volume of the positive electrode will involve one volume of the
% separator for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.

% Last CV in the positive electrode: derivative approximation with a central scheme
dCe_last_p = 2*(ce(param.Np+1)-ce(param.Np-1))/(3 * param.deltax_p*param.len_p + param.deltax_s*param.len_s);

% Separator

% Attention! The first volume of the separator will involve one volume of the
% positive section for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.

% First CV in the separator: derivative approximation with a central scheme
dCe_first_s = 2*(ce(param.Np+2)-ce(param.Np))/( param.deltax_p*param.len_p + 3* param.deltax_s*param.len_s);

% Central differentiation scheme
dCes = (ce(param.Np+3:param.Np+param.Ns)-ce(param.Np+1:param.Np+param.Ns-2))/(2*param.deltax_s*param.len_s);

% Attention! The last volume of the separator will involve one volume of the
% negative section for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.

% Last CV in the separator: derivative approximation with a central scheme
dCe_last_s = 2*(ce(param.Np+param.Ns+1)-ce(param.Np+param.Ns-1))/( param.deltax_n*param.len_n + 3*param.deltax_s*param.len_s);

% Negative electrode

% Attention! The first volume of the negative electrode will involve one volume of the
% separator section for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.

% First CV in the negative electrode: derivative approximation with a central scheme
dCe_first_n = 2*(ce(param.Np+param.Ns+2)-ce(param.Np+param.Ns))/(3 * param.deltax_n*param.len_n + param.deltax_s*param.len_s);

dCen = [(ce(param.Np+param.Ns+3:end)-ce(param.Np+param.Ns+1:end-2))/(2*param.deltax_p*param.len_p);... 	% Central differentiation scheme
    (3*ce(end)-4*ce(end-1)+ce(end-2))/(2*param.deltax_n*param.len_n) 							% Backward differentiation scheme
    ];

dCe = [dCep;dCe_last_p;dCe_first_s;dCes;dCe_last_s;dCe_first_n;dCen];

end
