%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLATEELECTROLYTECONDUCTIVITIES evaluates the interpolation of the
% electrolyte conductivities at the edges of the control volumes. The
% returned values are obtained using the harmonic mean.

function [Keff_p_medio, Keff_s_medio, Keff_n_medio] = interpolateElectrolyteConductivities(Keff_p,Keff_s,Keff_n,param)
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

Keff_n_medio = [Keff_n_medio;0];
end