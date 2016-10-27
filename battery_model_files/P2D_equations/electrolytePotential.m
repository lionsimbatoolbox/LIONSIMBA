%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELECTROLYTEPOTENTIAL evaluates the residuals for the electrolyte
% potential.

function [res_Phie,Keff] = electrolytePotential(jflux,ce,T,param,Phie)

%% Effective electrolyte conductivity
% Keff = param.ElectrolyteConductivityFunction(ce,T,param);
% keyboard
% Comment this for benchmark purposes
Keff_p = param.ElectrolyteConductivityFunction(ce(1:param.Np),T(param.Nal+1:param.Nal+param.Np),param,'p');
Keff_s = param.ElectrolyteConductivityFunction(ce(param.Np+1:param.Np+param.Ns),T(param.Nal+param.Np+1:param.Nal+param.Np+param.Ns),param,'s');
Keff_n = param.ElectrolyteConductivityFunction(ce(param.Np+param.Ns+1:end),T(param.Nal+param.Np+param.Ns+1:end-param.Nco),param,'n');

% Uncomment this for benchmark purposes
% Keff_p = param.eps_i(1)^param.brugg_p *(4.1253*1e-2 + 5.007*1e-4*ce(1:param.Np) - 4.7212*1e-7*ce(1:param.Np).^2 +1.5094*1e-10*ce(1:param.Np).^3 -1.6018*1e-14*ce(1:param.Np).^4);
% Keff_s = param.eps_i(2)^param.brugg_s *(4.1253*1e-2 + 5.007*1e-4*ce(param.Np+1:param.Np+param.Ns) - 4.7212*1e-7*ce(param.Np+1:param.Np+param.Ns).^2 +1.5094*1e-10*ce(param.Np+1:param.Np+param.Ns).^3 -1.6018*1e-14*ce(param.Np+1:param.Np+param.Ns).^4);
% Keff_n = param.eps_i(3)^param.brugg_n *(4.1253*1e-2 + 5.007*1e-4*ce(param.Np+param.Ns+1:end) - 4.7212*1e-7*ce(param.Np+param.Ns+1:end).^2 +1.5094*1e-10*ce(param.Np+param.Ns+1:end).^3 -1.6018*1e-14*ce(param.Np+param.Ns+1:end).^4);

        
Keff = [Keff_p;Keff_s;Keff_n];
[Keff_p_medio, Keff_s_medio, Keff_n_medio] = interpolateElectrolyteConductivities(Keff_p,Keff_s,Keff_n,param);
%% Matrix building

% i-th element
A_p = diag(Keff_p_medio);
A_p(2:end,2:end)    = A_p(2:end,2:end) + diag(Keff_p_medio(1:end-1));
% i+1-th element
A_p(1:end-1,2:end)  = A_p(1:end-1,2:end) - diag(Keff_p_medio(1:end-1));
% i-1-th element
A_p(2:end,1:end-1)  = A_p(2:end,1:end-1) - diag(Keff_p_medio(1:end-1));

% i-th element
A_s = diag(Keff_s_medio);
A_s(2:end,2:end)    = A_s(2:end,2:end) + diag(Keff_s_medio(1:end-1));
% i+1-th element
A_s(1:end-1,2:end)  = A_s(1:end-1,2:end) - diag(Keff_s_medio(1:end-1));
% i-1-th element
A_s(2:end,1:end-1)  = A_s(2:end,1:end-1) - diag(Keff_s_medio(1:end-1));

% i-th element
A_n = diag(Keff_n_medio);
A_n(2:end,2:end)    = A_n(2:end,2:end) + diag(Keff_n_medio(1:end-1));
% i+1-th element
A_n(1:end-1,2:end)  = A_n(1:end-1,2:end) - diag(Keff_n_medio(1:end-1));
% i-1-th element
A_n(2:end,1:end-1)  = A_n(2:end,1:end-1) - diag(Keff_n_medio(1:end-1));

A_n = A_n./(param.deltax_n*param.len_n);
A_s = A_s./(param.deltax_s*param.len_s);
A_p = A_p./(param.deltax_p*param.len_p);

A_tot = blkdiag(A_p,A_s,A_n);

% Fix values to enforce BC on the left side of the positive electrode.
A_tot(1,1:2) = [Keff_p_medio(1) -Keff_p_medio(1)]./(param.deltax_p*param.len_p);

% The value of Phie in the last volume of the negative electrode is known
% and fixed.
A_tot(end,end-1:end) = [0 1];

%% Interfaces Positive electrode (last volume of the positive)

% Here we are in the last volume of the positive
den_sn = (param.deltax_p*param.len_p/2+param.deltax_s*param.len_s/2);
last_p = Keff_p_medio(end-1)/(param.deltax_p*param.len_p);
A_tot(param.Np,param.Np-1:param.Np+1) = [-last_p (last_p+Keff_p_medio(end)/den_sn) -Keff_p_medio(end)/den_sn];
%% Interfaces Positive electrode (first volume of the separator)

% Here we are in the first volume of the separator
den_sp = (param.deltax_p*param.len_p/2+param.deltax_s*param.len_s/2);
first_s = Keff_s_medio(1)/(param.deltax_s*param.len_s);
A_tot(param.Np+1,param.Np:param.Np+2) = [-Keff_p_medio(end)/den_sp (first_s+Keff_p_medio(end)/den_sp) -first_s];

%% Interfaces Positive electrode (last volume of the separator)
% Here we are in the last volume of the separator
den_sn = (param.deltax_n*param.len_n/2+param.deltax_s*param.len_s/2);
last_s = Keff_s_medio(end-1)/(param.deltax_s*param.len_s);
A_tot(param.Np+param.Ns,param.Np+param.Ns-1:param.Np+param.Ns+1) = [-last_s (last_s+Keff_s_medio(end)/den_sn) -Keff_s_medio(end)/den_sn];

%% Interfaces Positive electrode (first volume of the negative)
% Here we are inside the first volume of the negative electrode
den_ns = (param.deltax_n*param.len_n/2+param.deltax_s*param.len_s/2);
first_n = Keff_n_medio(2)/(param.deltax_n*param.len_n);
A_tot(param.Np+param.Ns+1,param.Np+param.Ns:param.Np+param.Ns+2) = [-Keff_s_medio(end)/den_ns (first_n+Keff_s_medio(end)/den_sn) -first_n];
% Fix this value otherwise it will be not correct
A_tot(param.Np+param.Ns+2,param.Np+param.Ns+1) = -Keff_n_medio(1)/(param.deltax_n*param.len_n);

%% Electrolyte concentration interpolation
% Evaluate the interpolation of the electrolyte concentration values at the
% edges of the control volumes.
[ce_mean_p, ce_mean_ps, ce_mean_s, ce_mean_sn, ce_mean_n] = interpolateElectrolyteConcentration(ce,param);
%% Temperature interpolation
% Evaluate the temperature value at the edges of the control volumes
[T_mean_p, T_mean_ps, T_mean_s, T_mean_sn, T_mean_n] = interpolateTemperature(T,param);
%% Electrolyte fluxes
% Evaluate the interpolation of the electrolyte concentration fluxes at the
% edges of the control volumes.
[ce_flux_p, ce_flux_ps, ce_flux_s, ce_flux_sn, ce_flux_n] = interpolateElectrolyteConcetrationFluxes(ce,param);
%% RHS arrays
K = 2*param.R*(1-param.tplus) / param.F;

prod_p = (Keff_p_medio.*[T_mean_p;T_mean_ps].*[ce_flux_p;ce_flux_ps].*[1./ce_mean_p;1/ce_mean_ps]);

prod_s = (Keff_s_medio.*[T_mean_s;T_mean_sn].*[ce_flux_s;ce_flux_sn].*[1./ce_mean_s;1/ce_mean_sn]);

% The last element of Keff_n_medio is not needed
prod_n = (Keff_n_medio(1:end-1).*T_mean_n.*ce_flux_n.*1./ce_mean_n);

if(isa(prod_p,'casadi.SX') || isa(prod_p,'casadi.MX'))
    prod_tot = [prod_p;prod_s;prod_n];
    prod_tot = prod_tot(2:end)-prod_tot(1:end-1);
else
    prod_tot = diff([prod_p;prod_s;prod_n]);
end
prod_tot = [prod_p(1);prod_tot];
flux_p = param.deltax_p*param.len_p*param.F*param.a_i(1)*jflux(1:param.Np);

flux_s = param.deltax_s*param.len_s*param.F*param.a_i(2);

flux_n = param.deltax_n*param.len_n*param.F*param.a_i(3)*jflux(param.Np+1:end-1);

flux_tot = [flux_p;flux_s*ones(param.Ns,1);flux_n];

f = flux_tot-K*prod_tot;
% Set the last element of Phie to 0
% f(end+1) = 0;
f=[f;0];

if(~isa(prod_p,'casadi.SX') && ~isa(prod_p,'casadi.MX'))
    A_tot = sparse(A_tot);
end

res_Phie = A_tot*Phie-f;
end