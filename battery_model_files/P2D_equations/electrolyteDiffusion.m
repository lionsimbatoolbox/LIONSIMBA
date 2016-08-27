%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELECTROLYTEDIFFUSION evaluates the residual for the electrolyte
% concentration of Li-ions in the electrolyte solution.

function resCe = electrolyteDiffusion(t,ce,dCe,jflux,T,param)

% Diffusion coefficients
% Comment this for benchmark purposes
Deff_p = param.ElectrolyteDiffusionFunction(ce(1:param.Np),T(param.Nal+1:param.Nal+param.Np),param,'p');
Deff_s = param.ElectrolyteDiffusionFunction(ce(param.Np+1:param.Np+param.Ns),T(param.Nal+param.Np+1:param.Nal+param.Np+param.Ns),param,'s');
Deff_n = param.ElectrolyteDiffusionFunction(ce(param.Np+param.Ns+1:end),T(param.Nal+param.Np+param.Ns+1:end-param.Nco),param,'n');

% Uncomment this for benchmark purposes
% Deff_p = repmat(param.Dp*param.eps_i(1)^param.brugg_p,param.Np,1);
% Deff_s = repmat(param.Ds*param.eps_i(2)^param.brugg_s,param.Ns,1);
% Deff_n = repmat(param.Dn*param.eps_i(3)^param.brugg_n,param.Nn,1);

% Interpolation of the diffusion coefficients
[Deff_p, Deff_s, Deff_n] = interpolateDiffusionCoefficients(Deff_p,Deff_s,Deff_n,param);

%% Positive A Matrix
A_p = - diag(Deff_p);

A_p(2:end,2:end)    = A_p(2:end,2:end) - diag(Deff_p(1:end-1));
A_p(1:end-1,2:end)  = A_p(1:end-1,2:end) + diag(Deff_p(1:end-1));
A_p(2:end,1:end-1)  = A_p(2:end,1:end-1) + diag(Deff_p(1:end-1));
%% Separator A Matrix
A_s = - diag(Deff_s);

A_s(2:end,2:end)    = A_s(2:end,2:end) - diag(Deff_s(1:end-1));
A_s(1:end-1,2:end)  = A_s(1:end-1,2:end) + diag(Deff_s(1:end-1));
A_s(2:end,1:end-1)  = A_s(2:end,1:end-1) + diag(Deff_s(1:end-1));
%% Negative A Matrix
A_n = - diag(Deff_n);

A_n(2:end,2:end)    = A_n(2:end,2:end) - diag(Deff_n(1:end-1));
A_n(1:end-1,2:end)  = A_n(1:end-1,2:end) + diag(Deff_n(1:end-1));
A_n(2:end,1:end-1)  = A_n(2:end,1:end-1) + diag(Deff_n(1:end-1));
% Fix the last elements of the A_n
A_n(end,end-1:end)  = [Deff_n(end-1) -Deff_n(end-1)];
%% A_tot matrix

A_tot                   = blkdiag(A_p,A_s,A_n);

% Divide by the deltax and the length of the positive electrode
A_tot(1:param.Np,1:param.Np)    = A_tot(1:param.Np,1:param.Np)/(param.deltax_p^2*param.len_p^2);
% Reset values on the lines for the interfaces conditions
A_tot(param.Np,:)           = 0;
A_tot(param.Np+1,:)         = 0;

% Divide by the deltax and the length of the separator
A_tot(param.Np+1:param.Np+param.Ns,param.Np+1:param.Np+param.Ns)    = A_tot(param.Np+1:param.Np+param.Ns,param.Np+1:param.Np+param.Ns)/(param.deltax_s^2 * param.len_s^2);
% Reset values on the lines for the interfaces conditions
A_tot(param.Np+param.Ns,:)      = 0;
A_tot(param.Np+param.Ns+1,:)    = 0;

% Divide by the deltax and the length of the negative electrode
A_tot(param.Np+param.Ns+1:end,param.Np+param.Ns+1:end)      = A_tot(param.Np+param.Ns+1:end,param.Np+param.Ns+1:end)/(param.deltax_n^2 * param.len_n^2);

%% Interface between separator and positive electrode (last volume in the positive electrode)

% Compute the common denominator at the interface
den_s   = (param.deltax_p*param.len_p/2 + param.deltax_s*param.len_s/2);
% Last diffusion coefficient of the positive electrode
last_p  = Deff_p(end-1)/(param.deltax_p*param.len_p);
% Diffusion coefficient on the interface
first_s = Deff_p(end)/den_s;
% Fix the values at the boundaries
A_tot(param.Np,param.Np-1:param.Np+1) = [last_p -(last_p+ first_s) first_s]/(param.deltax_p*param.len_p*param.eps_i(1));

%% Interface between separator and positive electrode (first volume in the separator)

% Compute the common denominator at the interface
den_s       = (param.deltax_p*param.len_p/2 + param.deltax_s*param.len_s/2);
% First diffusion coefficient in the separator
second_s    = Deff_s(1)/(param.deltax_s*param.len_s);
% Diffusion coefficient on the interface
first_s     = Deff_p(end)/den_s;

A_tot(param.Np+1,param.Np:param.Np+2) = [first_s -(first_s+second_s) second_s]/(param.deltax_s*param.len_s*param.eps_i(2));

%% Interface between separator and negative electrode (last volume in the separator)

% Compute the common denominator at the interface
den_s   = (param.deltax_s*param.len_s/2 + param.deltax_n*param.len_n/2);
% Last diffusion coefficient in the separator
last_s  = Deff_s(end-1)/(param.deltax_s*param.len_s);
% Diffusion coefficient on the interface
first_n = Deff_s(end)/den_s;

A_tot(param.Np+param.Ns,param.Np+param.Ns-1:param.Np+param.Ns+1) = [last_s -(last_s+first_n) first_n]/(param.deltax_s*param.len_s*param.eps_i(2));

%% Interface between separator and negative electrode (first volume in the negative electrode)

% Compute the common denominator at the interface
den_n       = (param.deltax_s*param.len_s/2 + param.deltax_n*param.len_n/2);
% First diffusion coefficient in the negative electrode
second_n    = Deff_n(1)/(param.deltax_n*param.len_n);
% Diffusion coefficient on the interface
first_n     = Deff_s(end)/den_n;

A_tot(param.Np+param.Ns+1,param.Np+param.Ns:param.Np+param.Ns+2) = [first_n -(first_n+second_n) second_n]/(param.deltax_n*param.len_n*param.eps_i(3));

%% Usefull stuff
a_tot       = [
    repmat(param.a_i(1),param.Np,1);...
    zeros(param.Ns,1);...
    repmat(param.a_i(3),param.Nn,1)...
    ];

jflux_tot   = [
    jflux(1:param.Np);...
    zeros(param.Ns,1);...
    jflux(param.Np+1:end)...
    ];

eps_tot     = [
    repmat(param.eps_i(1),param.Np,1);...
    repmat(param.eps_i(2),param.Ns,1);...
    repmat(param.eps_i(3),param.Nn,1)
    ];

% Build porosities matrix
K = 1./eps_tot;
A_eps = diag(K);
A_eps = A_eps + diag(K(1:end-1),1);
A_eps = A_eps + diag(K(1:end-1),-1);
A_eps = sparse(A_eps);

A_eps(param.Np,param.Np-1:param.Np+1) = 1;
A_eps(param.Np+1,param.Np:param.Np+2) = 1;

A_eps(param.Np+param.Ns,param.Np+param.Ns-1:param.Np+param.Ns+1) = 1;
A_eps(param.Np+param.Ns+1,param.Np+param.Ns:param.Np+param.Ns+2) = 1;

G = A_eps.*A_tot;

resCe = dCe - (G*ce + K.*(1-param.tplus).*a_tot.*jflux_tot);
end