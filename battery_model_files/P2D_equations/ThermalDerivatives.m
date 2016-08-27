%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THERMALDERIVATIVES evaluetes the set of derivatives used in the thermal
% dynamics.

function [dPhis, dPhie, dCe, ie, is, die, dis] = ThermalDerivatives(Phis,Phie,ce,param)
% This script evaluates the derivatives of the dependent variables for
% thermal calculations. Note that the derivatives are evaluated within the
% control volume and not at the boundaries.

%% Solid potential derivatives

% The derivatives in the first and last volumes are evaluated with first
% order approximations, while the mid points ones with second order
% approximations.

% Positive Electrode
dPhisp = [(-3*Phis(1)+4*Phis(2)-Phis(3))/(2*param.deltax_p*param.len_p);...
    (Phis(3:param.Np)-Phis(1:param.Np-2)) / (2*param.deltax_p*param.len_p);...
    (3*Phis(param.Np)-4*Phis(param.Np-1)+Phis(param.Np-2)) / (2*param.deltax_p*param.len_p);...
    ];

% Negative Electrode
dPhisn = [(-3*Phis(param.Np+1)+4*Phis(param.Np+2)-Phis(param.Np+3))/(2*param.deltax_n*param.len_n);...
    (Phis(param.Np+3:end)-Phis(param.Np+1:end-2)) / (2*param.deltax_n*param.len_n);...
    (3*Phis(end)-4*Phis(end-1)+Phis(end-2)) / (2*param.deltax_n*param.len_n);...
    ];

dPhis = [dPhisp;dPhisn];

%% Electrolyte potential

% Positive Electrode

dPhiep = [ (-3*Phie(1)+4*Phie(2)-Phie(3))/(2*param.deltax_p*param.len_p);...
    (Phie(3:param.Np)-Phie(1:param.Np-2))/(2*param.deltax_p*param.len_p);...
    ];

% Attention! The last volume of the positive will involve one volume of the
% separator for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.

% Last volume derivative
dPhie_last_p = (Phie(param.Np+1)-Phie(param.Np-1))/(3 * param.deltax_p*param.len_p + param.deltax_s*param.len_s);

% Separator

% Attention! The first volume of the separator will involve one volume of the
% positive section for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.

dPhie_first_s = (Phie(param.Np+2)-Phie(param.Np))/(param.deltax_p*param.len_p + 3* param.deltax_s*param.len_s);

dPhies =  (Phie(param.Np+3:param.Np+param.Ns)-Phie(param.Np+1:param.Np+param.Ns-2))/(2*param.deltax_s*param.len_s);

% Attention! The last volume of the separator will involve one volume of the
% negative section for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.

dPhie_last_s = (Phie(param.Np+param.Ns+1)-Phie(param.Np+param.Ns-1))/( param.deltax_n*param.len_n + 3*param.deltax_s*param.len_s);

% Negative electrode

% Attention! The first volume of the negative will involve one volume of the
% separator section for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.
dPhie_first_n = (Phie(param.Np+param.Ns+2)-Phie(param.Np+param.Ns))/(3 * param.deltax_n*param.len_n + param.deltax_s*param.len_s);

dPhien = [(Phie(param.Np+param.Ns+3:end)-Phie(param.Np+param.Ns+1:end-2))/(2*param.deltax_p*param.len_p);...
            (3*Phie(end)-4*Phie(end-1)+Phie(end-2))/(2*param.deltax_n*param.len_n)
    ];

dPhie = [dPhiep;dPhie_last_p;dPhie_first_s;dPhies;dPhie_last_s;dPhie_first_n;dPhien];

%% Electrolyte concentration

% Positive Electrode

dCep = [ (-3*ce(1)+4*ce(2)-ce(3))/(2*param.deltax_p*param.len_p);...
    (ce(3:param.Np)-ce(1:param.Np-2))/(2*param.deltax_p*param.len_p);...
    ];

% Attention! The last volume of the positive will involve one volume of the
% separator for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.

% Last volume derivative
dCe_last_p = (ce(param.Np+1)-ce(param.Np-1))/(3 * param.deltax_p*param.len_p + param.deltax_s*param.len_s);

% Separator

% Attention! The first volume of the separator will involve one volume of the
% positive section for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.

dCe_first_s = (ce(param.Np+2)-ce(param.Np))/( param.deltax_p*param.len_p + 3* param.deltax_s*param.len_s);

dCes = (ce(param.Np+3:param.Np+param.Ns)-ce(param.Np+1:param.Np+param.Ns-2))/(2*param.deltax_s*param.len_s);

% Attention! The last volume of the separator will involve one volume of the
% negative section for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.

dCe_last_s = (ce(param.Np+param.Ns+1)-ce(param.Np+param.Ns-1))/( param.deltax_n*param.len_n + 3*param.deltax_s*param.len_s);

% Negative electrode

% Attention! The first volume of the negative will involve one volume of the
% separator section for the calculation of the derivative. Therefore suitable
% considerations must be done with respect to the deltax quantities.
dCe_first_n = (ce(param.Np+param.Ns+2)-ce(param.Np+param.Ns))/(3 * param.deltax_n*param.len_n + param.deltax_s*param.len_s);

dCen = [(ce(param.Np+param.Ns+3:end)-ce(param.Np+param.Ns+1:end-2))/(2*param.deltax_p*param.len_p);...
            (3*ce(end)-4*ce(end-1)+ce(end-2))/(2*param.deltax_n*param.len_n)
    ];

dCe = [dCep;dCe_last_p;dCe_first_s;dCes;dCe_last_s;dCe_first_n;dCen];


%% Currents

iep = dPhis(1:param.Np)*param.sig_eff(1)+param.I;
isp = param.I-iep;
ien = dPhis(param.Np+1:end)*param.sig_eff(3)+param.I;
isn = param.I-ien;

diep = [(-3*iep(1)+4*iep(2)-iep(3))/(2*param.len_p*param.deltax_p);...
    (iep(3:end)-iep(1:end-2))/(2*param.len_p*param.deltax_p);...
    (3*iep(end)-4*iep(end-1)+iep(end-2))/(2*param.len_p*param.deltax_p)];

disp = [(-3*isp(1)+4*isp(2)-isp(3))/(2*param.len_p*param.deltax_p);...
    (isp(3:end)-isp(1:end-2))/(2*param.len_p*param.deltax_p);...
    (3*isp(end)-4*isp(end-1)+isp(end-2))/(2*param.len_p*param.deltax_p)];

dien = [(-3*ien(1)+4*ien(2)-ien(3))/(2*param.len_n*param.deltax_n);...
    (ien(3:end)-ien(1:end-2))/(2*param.len_n*param.deltax_n);...
    (+3*ien(end)-4*ien(end-1)+ien(end-2))/(2*param.len_n*param.deltax_n)];

disn = [(-3*isn(1)+4*isn(2)-isn(3))/(2*param.len_n*param.deltax_n);...
    (isn(3:end)-isn(1:end-2))/(2*param.len_n*param.deltax_n);...
    (+3*isn(end)-4*isn(end-1)+isn(end-2))/(2*param.len_n*param.deltax_n)];

ie=[iep;ien];
is=[isp;isn];
die=[diep;dien];
dis=[disp;disn];

end