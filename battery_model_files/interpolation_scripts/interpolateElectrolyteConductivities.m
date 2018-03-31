function [Keff_p_medio, Keff_s_medio, Keff_n_medio] = interpolateElectrolyteConductivities(Keff_p,Keff_s,Keff_n,param)
%	interpolateElectrolyteConductivities interpolates electrolyte conductivities at the edges of control volumes using harmonic mean.

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

%% Positive electrode mean conductivity
beta_p = 0.5;
Keff_p_medio = Keff_p(1:end-1).*Keff_p(2:end)./(beta_p*Keff_p(2:end)+(1-beta_p)*Keff_p(1:end-1));

% The last element of Keff_p_medio will be the harmonic mean of the
% elements at the interface positive-separator

beta_p_s = param.deltax_p*param.len_p/2 /(param.deltax_s*param.len_s/2+param.deltax_p*param.len_p/2);

Keff_p_s_interface = Keff_p(end)*Keff_s(1) / (beta_p_s*Keff_s(1) + (1-beta_p_s)*Keff_p(end));

Keff_p_medio = [Keff_p_medio;Keff_p_s_interface];

%% Separator mean conductivity
% Compute the harmonic mean values for the separator
beta_s = 0.5;
Keff_s_medio = Keff_s(1:end-1).*Keff_s(2:end)./(beta_s*Keff_s(2:end)+(1-beta_s)*Keff_s(1:end-1));

% The last element of Keff_s_medio will be the harmonic mean of the
% elements at the interface separator-negative

beta_s_n = param.deltax_s*param.len_s/2 /(param.deltax_s*param.len_s/2+param.deltax_n*param.len_n/2);

Keff_s_n_interface = Keff_s(end)*Keff_n(1) / (beta_s_n*Keff_n(1) + (1-beta_s_n)*Keff_s(end));

Keff_s_medio = [Keff_s_medio;Keff_s_n_interface];

%% Negative electrode mean conductivity
% Compute the harmonic mean values for the negative electrode

beta_n = 0.5;
Keff_n_medio = Keff_n(1:end-1).*Keff_n(2:end)./(beta_n*Keff_n(2:end)+(1-beta_n)*Keff_n(1:end-1));

Keff_n_medio = [Keff_n_medio;0]; % the cc interface is not used. The zero is only to match the dimensions
end