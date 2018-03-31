function [Deff_p_medio, Deff_s_medio, Deff_n_medio] = interpolateDiffusionCoefficients(Deff_p,Deff_s,Deff_n,param)
%	interpolateDiffusionCoefficients Evaluates the interpolation of the
%	values for the diffusion coefficients.

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

%% Diffusion coefficients for the positive electrode

% Beta coefficient for the interface with the separator
betaD_ps        = param.deltax_p*param.len_p/2 / (param.deltax_p*param.len_p/2 + param.deltax_s*param.len_s/2);
% Beta coefficient within the positive electrode
betaD_p        = 0.5;
% Harmonic mean for the coefficients within the positive electrode
Deff_p_medio    = (Deff_p(1:end-1).*Deff_p(2:end)./(betaD_p*Deff_p(2:end)+(1-betaD_p)*Deff_p(1:end-1)));
% Harmonic mean for the coefficients at the interface with the separator
Deff_medio_ps   = Deff_p(end)*Deff_s(1) / (betaD_ps*Deff_s(1)+(1-betaD_ps)*Deff_p(end));
% Build the overall array
Deff_p_medio    = [Deff_p_medio;Deff_medio_ps];

%%  Diffusion coefficients for the separator
% Beta coefficient for the interface with the negative electrode
betaD_sn = param.deltax_s*param.len_s/2 / (param.deltax_n*param.len_n/2 + param.deltax_s*param.len_s/2);
% Beta coefficient within the separator
betaD_s        = 0.5;
% Harmonic mean for the coefficients within the separator
Deff_s_medio = (Deff_s(1:end-1).*Deff_s(2:end))./(betaD_s*Deff_s(2:end)+(1-betaD_s)*Deff_s(1:end-1));
% Harmonic mean for the coefficients at the interface with the negative
% electrode
Deff_medio_sn = Deff_n(1)*Deff_s(end) / (betaD_sn*Deff_n(1)+(1-betaD_sn)*Deff_s(end));
% Build the overall array
Deff_s_medio = [Deff_s_medio;Deff_medio_sn];

%% Diffusion coefficients for the negative electrode
% Beta coefficient within the separator
betaD_n        = 0.5;
% Build the overall array. (Note that a 0 is added at the end of the array
% in order to have the correct dimension; in the end the 0 will be
% replaced with the right values in the A_tot matrix).
Deff_n_medio = [(Deff_n(1:end-1).*Deff_n(2:end))./(betaD_n*Deff_n(2:end)+(1-betaD_n)*Deff_n(1:end-1));0];

end
