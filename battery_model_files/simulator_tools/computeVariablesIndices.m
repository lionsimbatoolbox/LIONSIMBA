%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function param = computeVariablesIndices(param)
% Store the indices of the differential variables
param.ce_indices         = (1:param.Nsum);

% Modify the solid phase indices according to the model used. If Fick's
% law is used, then it is necessary to account also for the diffusion
% inside the solid particles.
if(param.SolidPhaseDiffusion==1 || param.SolidPhaseDiffusion==2)
    param.cs_average_indices    = (param.ce_indices(end)+1:param.ce_indices(end)+param.Np+param.Nn);
elseif(param.SolidPhaseDiffusion==3)
    param.cs_average_indices = (param.ce_indices(end)+1:param.ce_indices(end)+param.Np*param.Nr_p+param.Nn*param.Nr_n);
end

param.T_indices          = (param.cs_average_indices(end)+1:param.cs_average_indices(end)+param.Nal+param.Nsum+param.Nco);
param.film_indices       = (param.T_indices(end)+1:param.T_indices(end)+param.Nn);
param.Q_indices          = (param.film_indices(end)+1:param.film_indices(end)+param.Np+param.Nn);

% Store the indices of the algebraic variables.
param.jflux_indices      = (param.Q_indices(end)+1:param.Q_indices(end)+param.Np+param.Nn);
param.Phis_indices       = (param.jflux_indices(end)+1:param.jflux_indices(end)+param.Np+param.Nn);
param.Phie_indices       = (param.Phis_indices(end)+1:param.Phis_indices(end)+param.Np+param.Ns+param.Nn);
param.js_indices         = (param.Phie_indices(end)+1:param.Phie_indices(end)+param.Nn);
param.Iapp_indices       = (param.js_indices(end)+1);

end