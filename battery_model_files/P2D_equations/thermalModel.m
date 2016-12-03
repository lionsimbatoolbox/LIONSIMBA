%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THERMALMODEL evaluates the set of equations for the thermal dynamics.
function [res_dT, rhsT] = thermalModel(ce,Phie,Phis,Keff,jflux,T,dT,Up,Un,dudt_p,dudt_n,param)

% Effective electrolyte conductivities
Keff_p = Keff(1:param.Np);
Keff_s = Keff(param.Np+1:param.Np+param.Ns);
Keff_n = Keff(param.Np+param.Ns+1:end);

% Evaluate the derivatives used in Qohm calculations
[dPhis, dPhie, dCe] = ThermalDerivatives(Phis,Phie,ce,param);

%% Reversible heat generation rate

% Positive electrode
Qrev_p = param.F*param.a_i(1)*jflux(1:param.Np).*T(param.Nal+1:param.Nal+param.Np).*dudt_p;

% Negative Electrode
Qrev_n = param.F*param.a_i(3)*jflux(param.Np+1:end).*T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn).*dudt_n;

%% Reaction heat generation rate

% Positive overpotential
eta_p = (Phis(1:param.Np)-Phie(1:param.Np)-Up);
% Positive reaction heat generation rate
Qrxn_p = param.F*param.a_i(1)*jflux(1:param.Np).*eta_p;

% Negative overpotential
eta_n = (Phis(param.Np+1:end)-Phie(param.Np+param.Ns+1:end)-Un);
% Negative reaction heat generation rate
Qrxn_n = param.F*param.a_i(3)*jflux(param.Np+1:end).*eta_n;

%% Ohmic heat generation rate

% Positive electrode ohmic generation rate
Qohm_p = param.sig_eff(1) * (dPhis(1:param.Np)).^2 + Keff_p.*(dPhie(1:param.Np)).^2 + 2*param.R*Keff_p.*T(param.Nal+1:param.Nal+param.Np)*(1-param.tplus)/param.F.*dCe(1:param.Np).*1./ce(1:param.Np).*dPhie(1:param.Np);
% Separator ohmic generation rate
Qohm_s = Keff_s.*(dPhie(param.Np+1:param.Np+param.Ns)).^2 + 2*param.R*Keff_s.*T(param.Nal+param.Np+1:param.Nal+param.Np+param.Ns)*(1-param.tplus)/param.F.*dCe(param.Np+1:param.Np+param.Ns).*1./ce(param.Np+1:param.Np+param.Ns).*dPhie(param.Np+1:param.Np+param.Ns);
% Negative electrode ohmic generation rate
Qohm_n = param.sig_eff(3) * (dPhis(param.Np+1:end)).^2 +Keff_n.*(dPhie(param.Np+param.Ns+1:end)).^2 + 2*param.R*Keff_n.*T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn)*(1-param.tplus)/param.F.*dCe(param.Np+param.Ns+1:end).*1./ce(param.Np+param.Ns+1:end).*dPhie(param.Np+param.Ns+1:end);

%% A Matrix definitions

% Positive current collector
Lambda_al   = param.Lambda_al * ones(param.Nal,1);
A_al        = - diag(Lambda_al);

A_al(2:end,2:end)    = A_al(2:end,2:end) - diag(Lambda_al(1:end-1));
A_al(1:end-1,2:end)  = A_al(1:end-1,2:end) + diag(Lambda_al(1:end-1));
A_al(2:end,1:end-1)  = A_al(2:end,1:end-1) + diag(Lambda_al(1:end-1));

% Positive electrode
Lambda_p    = param.Lambda_p * ones(param.Np,1);
A_p         = - diag(Lambda_p);

A_p(2:end,2:end)    = A_p(2:end,2:end) - diag(Lambda_p(1:end-1));
A_p(1:end-1,2:end)  = A_p(1:end-1,2:end) + diag(Lambda_p(1:end-1));
A_p(2:end,1:end-1)  = A_p(2:end,1:end-1) + diag(Lambda_p(1:end-1));

% Separator
Lambda_s    = param.Lambda_s * ones(param.Ns,1);
A_s         = - diag(Lambda_s);

A_s(2:end,2:end)    = A_s(2:end,2:end) - diag(Lambda_s(1:end-1));
A_s(1:end-1,2:end)  = A_s(1:end-1,2:end) + diag(Lambda_s(1:end-1));
A_s(2:end,1:end-1)  = A_s(2:end,1:end-1) + diag(Lambda_s(1:end-1));

% Negative electrode
Lambda_n    = param.Lambda_n * ones(param.Nn,1);
A_n         = - diag(Lambda_n);

A_n(2:end,2:end)    = A_n(2:end,2:end) - diag(Lambda_n(1:end-1));
A_n(1:end-1,2:end)  = A_n(1:end-1,2:end) + diag(Lambda_n(1:end-1));
A_n(2:end,1:end-1)  = A_n(2:end,1:end-1) + diag(Lambda_n(1:end-1));


% Negative current collector
Lambda_co   = param.Lambda_co * ones(param.Nco,1);
A_co        = - diag(Lambda_co);

A_co(2:end,2:end)    = A_co(2:end,2:end) - diag(Lambda_co(1:end-1));
A_co(1:end-1,2:end)  = A_co(1:end-1,2:end) + diag(Lambda_co(1:end-1));
A_co(2:end,1:end-1)  = A_co(2:end,1:end-1) + diag(Lambda_co(1:end-1));

A_co(end,end-1:end) = [param.Lambda_co -param.Lambda_co];

% Divide the matrices for their deltax

A_al    = A_al/(param.deltax_al*param.len_al)^2;
A_p     = A_p/(param.deltax_p*param.len_p)^2;
A_s     = A_s/(param.deltax_s*param.len_s)^2;
A_n     = A_n/(param.deltax_n*param.len_n)^2;
A_co    = A_co/(param.deltax_co*param.len_co)^2;

A_tot = blkdiag(A_al,A_p,A_s,A_n,A_co);

%% Interfaces

% Interface between aluminium current collector and positive electrode. We
% are in the last volume of the current collector
beta_al_p   = (param.deltax_al*param.len_al/2)/(param.deltax_al*param.len_al/2+param.deltax_p*param.len_p/2);
Lambda_al_p = param.Lambda_al * param.Lambda_p /(beta_al_p*param.Lambda_p + (1-beta_al_p)*param.Lambda_al);
den_al_p    = param.deltax_p*param.len_p/2 +param.deltax_al*param.len_al/2;
last_al     = param.Lambda_al / (param.deltax_al*param.len_al);
first_p     = Lambda_al_p/den_al_p;

A_tot(param.Nal,param.Nal-1:param.Nal+1) = [last_al -(last_al+first_p) first_p]/(param.deltax_al*param.len_al);

% Interface between aluminium current collector and positive electrode. We
% are in the first volume of the positive electrode

den_al_p    = param.deltax_p*param.len_p/2 +param.deltax_al*param.len_al/2;
second_p    = param.Lambda_p / (param.deltax_p*param.len_p);
first_p     = Lambda_al_p/den_al_p;

A_tot(param.Nal+1,param.Nal:param.Nal+2) = [first_p -(second_p+first_p) second_p]/(param.deltax_p*param.len_p);

% Interface between positive electrode and separator. We
% are in the last volume of the positive electrode
beta_p_s    = (param.deltax_p*param.len_p/2)/(param.deltax_s*param.len_s/2+param.deltax_p*param.len_p/2);
Lambda_p_s  = param.Lambda_s * param.Lambda_p /(beta_p_s*param.Lambda_s + (1-beta_p_s)*param.Lambda_p);

den_p_s     = param.deltax_p*param.len_p/2 +param.deltax_s*param.len_s/2;
last_p      = param.Lambda_p / (param.deltax_p*param.len_p);
first_s     = Lambda_p_s/den_p_s;

A_tot(param.Nal+param.Np,param.Nal+param.Np-1:param.Nal+param.Np+1) = [last_p -(last_p+first_s) first_s]/(param.deltax_p*param.len_p);


% Interface between positive electrode and separator. We
% are in the first volume of the separator
den_p_s     = param.deltax_p*param.len_p/2 +param.deltax_s*param.len_s/2;
second_s    = param.Lambda_s / (param.deltax_s*param.len_s);
first_s     = Lambda_p_s/den_p_s;

A_tot(param.Nal+param.Np+1,param.Nal+param.Np:param.Nal+param.Np+2) = [first_s -(second_s+first_s) second_s]/(param.deltax_s*param.len_s);

% Interface between separator negative electrode. We
% are in the last volume of the separator
beta_s_n    = (param.deltax_s*param.len_s/2)/(param.deltax_s*param.len_s/2+param.deltax_n*param.len_n/2);
Lambda_s_n  = param.Lambda_s * param.Lambda_n /(beta_s_n*param.Lambda_n + (1-beta_s_n)*param.Lambda_s);

den_s_n     = param.deltax_n*param.len_n/2 +param.deltax_s*param.len_s/2;
last_s      = param.Lambda_s / (param.deltax_s*param.len_s);
first_n     = Lambda_s_n/den_s_n;

A_tot(param.Nal+param.Np+param.Ns,param.Nal+param.Np+param.Ns-1:param.Nal+param.Np+param.Ns+1) = [last_s -(last_s+first_n) first_n]/(param.deltax_s*param.len_s);

% Interface between separator negative electrode. We
% are in the first volume of the negative electrode

den_s_n     = param.deltax_n*param.len_n/2 +param.deltax_s*param.len_s/2;
second_n    = param.Lambda_n / (param.deltax_n*param.len_n);
first_n     = Lambda_s_n/den_s_n;

A_tot(param.Nal+param.Np+param.Ns+1,param.Nal+param.Np+param.Ns:param.Nal+param.Np+param.Ns+2) = [first_n -(first_n+second_n) second_n]/(param.deltax_n*param.len_n);


% Interface between negative electrode and negative current collector. We
% are in the last volume of the negative electrode
beta_n_co   = (param.deltax_n*param.len_n/2)/(param.deltax_co*param.len_co/2+param.deltax_n*param.len_n/2);
Lambda_n_co = param.Lambda_co * param.Lambda_n /(beta_n_co*param.Lambda_co + (1-beta_n_co)*param.Lambda_n);

den_n_co    = param.deltax_n*param.len_n/2 +param.deltax_co*param.len_co/2;
last_n      = param.Lambda_n / (param.deltax_n*param.len_n);
first_co    = Lambda_n_co/den_n_co;

A_tot(param.Nal+param.Np+param.Ns+param.Nn,param.Nal+param.Np+param.Ns+param.Nn-1:param.Nal+param.Np+param.Ns+param.Nn+1) = [last_n -(last_n+first_co) first_co]/(param.deltax_n*param.len_n);


% Interface between negative electrode and negative current collector. We
% are in the last volume of the negative electrode

den_n_co    = param.deltax_n*param.len_n/2 +param.deltax_co*param.len_co/2;
second_co   = param.Lambda_co / (param.deltax_co*param.len_co);
first_co    = Lambda_n_co/den_n_co;

A_tot(param.Nal+param.Np+param.Ns+param.Nn+1,param.Nal+param.Np+param.Ns+param.Nn:param.Nal+param.Np+param.Ns+param.Nn+2) = [-first_co (second_co+first_co) -second_co]/(param.deltax_co*param.len_co);

if(~isa(T,'casadi.SX') && ~isa(T,'casadi.MX'))
    A_tot = sparse(A_tot);
    Qrev_tot = (sparse([zeros(param.Nal,1);Qrev_p;zeros(param.Ns,1);Qrev_n;zeros(param.Nco,1)]));
    Qrxn_tot = (sparse([zeros(param.Nal,1);Qrxn_p;zeros(param.Ns,1);Qrxn_n;zeros(param.Nco,1)]));
else
    Qrev_tot = (([zeros(param.Nal,1);Qrev_p;zeros(param.Ns,1);Qrev_n;zeros(param.Nco,1)]));
    Qrxn_tot = (([zeros(param.Nal,1);Qrxn_p;zeros(param.Ns,1);Qrxn_n;zeros(param.Nco,1)]));
end

Qohm_tot = [param.I^2/param.sig_al*ones(param.Nal,1);Qohm_p;Qohm_s;Qohm_n;param.I^2/param.sig_co*ones(param.Nco,1)];

BC = [param.T_BC_sx;zeros(param.Nal-1+param.Np+param.Ns+param.Nn+param.Nco-1,1);param.T_BC_dx];

predT = [param.rho_al*param.Cpal*ones(param.Nal,1);...
    param.rho_p*param.Cpp*ones(param.Np,1);...
    param.rho_s*param.Cps*ones(param.Ns,1);...
    param.rho_n*param.Cpn*ones(param.Nn,1);...
    param.rho_co*param.Cpco*ones(param.Nco,1);...
    ];

rhsT    = (1./predT.*(A_tot*T) + 1./predT.*Qrev_tot + 1./predT.*Qrxn_tot + 1./predT.*Qohm_tot + 1./predT.*BC);
res_dT  = dT - rhsT;


end