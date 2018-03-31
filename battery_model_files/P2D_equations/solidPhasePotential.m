function [res_Phis] = solidPhasePotential(jflux,param,Phis,I_density)
% solidPhasePotential evaluates the residuals of the solid potential equation.

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

%% Positive electrode

% RHS for the solid potential in the positive electrode. The BC on the left is enforced (Neumann BC)
f_p = ((param.len_p*param.deltax_p*param.a_i(1)*param.F*jflux(1))-I_density)*param.deltax_p*param.len_p/param.sig_eff(1);
% RHS for the solid potential in the positive electrode.
f_p = [f_p;(param.len_p^2*param.deltax_p^2*param.a_i(1)*param.F*jflux(2:param.Np))/param.sig_eff(1)];

%% Negative electrode

% RHS for the solid potential in the negative electrode.
f_n = (param.len_n^2*param.deltax_n^2*param.a_i(3)*param.F*jflux(param.Np+1:end-1))/param.sig_eff(3);

if param.OperatingMode==1 || param.OperatingMode==4 || param.OperatingMode==3 %The Neumann BC on the right is enforced only when operating using applied current density as the input
    % RHS for the solid potential in the negative electrode.
    f_n =[f_n;((param.len_n*param.deltax_n*param.a_i(3)*param.F*jflux(end))+I_density)*param.deltax_n*param.len_n/param.sig_eff(3)];
end

%% Residual array
% Return the residual array
res_Phis = [
    param.A_p*Phis(1:param.Np)-f_p;...
    param.A_n*Phis(param.Np+1:end)-f_n;
    ];
end
