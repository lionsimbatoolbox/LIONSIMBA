%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLIDPHASEPOTENTIAL evaluates the residuals of the solid potential equation. 

function [res_Phis] = solidPhasePotential(jflux,param,Phis)

%% Positive electrode


% RHS for the solid potential in the positive electrode. The BC on the left
% is enforced
f_p(1,1) = ((param.len_p*param.deltax_p*param.a_i(1)*param.F*jflux(1))-param.I)*param.deltax_p*param.len_p/param.sig_eff(1);
% RHS for the solid potential in the positive electrode.
f_p = [f_p;(param.len_p^2*param.deltax_p^2*param.a_i(1)*param.F*jflux(2:param.Np))/param.sig_eff(1)];

%% Negative electrode


% RHS for the solid potential in the negative electrode.
f_n = (param.len_n^2*param.deltax_n^2*param.a_i(3)*param.F*jflux(param.Np+1:end-1))/param.sig_eff(3);
% RHS for the solid potential in the negative electrode. The BC on the
% right is enforced
f_n =[f_n;((param.len_n*param.deltax_n*param.a_i(3)*param.F*jflux(end))+param.I)*param.deltax_n*param.len_n/param.sig_eff(3)];

%% Residual array
% Return the residual array
res_Phis = [
    param.A_p*Phis(1:param.Np)-f_p;...
    param.A_n*Phis(param.Np+1:end)-f_n;
    ];
end