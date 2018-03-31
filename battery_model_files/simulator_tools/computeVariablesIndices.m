function param = computeVariablesIndices(param)
%   computeVariablesIndices computes the indices of the differential variables
%   used in the code

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

param.ce_indices         = (1:param.Nsum);

% Modify the solid phase indices according to the model used. If Fick's
% law is used, then it is necessary to account also for the diffusion
% inside the solid particles.
if(param.SolidPhaseDiffusion==1 || param.SolidPhaseDiffusion==2)
    param.cs_average_indices    = (param.ce_indices(end)+1:param.ce_indices(end)+param.Np+param.Nn);
elseif(param.SolidPhaseDiffusion==3)
    param.cs_average_indices = (param.ce_indices(end)+1:param.ce_indices(end)+param.Np*param.Nr_p+param.Nn*param.Nr_n);
end

param.T_indices    = (param.cs_average_indices(end)+1:param.cs_average_indices(end)+param.Nal+param.Nsum+param.Ncu);
param.film_indices = (param.T_indices(end)+1:param.T_indices(end)+param.Nn);
param.Q_indices    = (param.film_indices(end)+1:param.film_indices(end)+param.Np+param.Nn);

% Store the indices of the algebraic variables.
param.jflux_indices     = (param.Q_indices(end)+1:param.Q_indices(end)+param.Np+param.Nn);
param.Phis_indices      = (param.jflux_indices(end)+1:param.jflux_indices(end)+param.Np+param.Nn);
param.Phie_indices      = (param.Phis_indices(end)+1:param.Phis_indices(end)+param.Np+param.Ns+param.Nn);
param.js_indices        = (param.Phie_indices(end)+1:param.Phie_indices(end)+param.Nn);
param.curr_dens_indices = (param.js_indices(end)+1);

end
