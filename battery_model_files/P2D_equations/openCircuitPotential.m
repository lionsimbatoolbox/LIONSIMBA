%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPENCIRCUITPOTENTIAL  Evaluates the Open Circuit Voltage of the cathode 
%                       and anode of a particular battery. 
%                       The measurement unit is [V]
%
%   [U_p,dudt_p,U_n,dudt_n] = OPENCIRCUITPOTENTIAL(cs_star,T,param)
%
%       - output:
%               - U_p and U_n       :  represents respectively the OCV of the 
%                                      anode and the cathode.
%               - dudt_p and dudt_n :  represents respectively the OCV 
%                                      variation with respect to temperature variations. 
%                                      This is considered only in temperature enabled simulations.
%
%
%   You can modify the way that the OCV are computed, as
%   function of electrolyte concentration and temperature. The main script
%   will pass also the param array.

function [U_p,dudt_p,U_n,dudt_n] = openCircuitPotential(cs_star,T,param)
	% Calculate the open circuit voltage of the battery in the positive
	% electrode
	theta_p  = cs_star(1:param.Np)./param.cs_maxp;

	% Compute the variation of OCV with respect to temperature variations [V/K]
	dudt_p   = -0.001 * (0.199521039-0.928373822*theta_p + 1.364550689000003*theta_p.^2-0.6115448939999998*theta_p.^3);
	dudt_p   = dudt_p./(1-5.661479886999997*theta_p +11.47636191*theta_p.^2-9.82431213599998*theta_p.^3+3.048755063*theta_p.^4);

	% Define the OCV for the positive electrode
	U_p    = (-4.656+88.669*theta_p.^2 - 401.119*theta_p.^4 + 342.909*theta_p.^6 - 462.471*theta_p.^8 + 433.434*theta_p.^10);
	U_p      = U_p./(-1+18.933*theta_p.^2-79.532*theta_p.^4+37.311*theta_p.^6-73.083*theta_p.^8+95.96*theta_p.^10);
	U_p      = U_p + (T(param.Nal+1:param.Nal+param.Np)-param.Tref).*dudt_p;


	% Calculate the open circuit voltage of the battery in the negative
	% electrode
	theta_n  = cs_star(param.Np+1:end)./ param.cs_maxn;

	% Compute the variation of OCV with respect to temperature variations [V/K]
	dudt_n = 0.001*(0.005269056 +3.299265709*theta_n-91.79325798*theta_n.^2+1004.911008*theta_n.^3-5812.278127*theta_n.^4 + ...
		19329.7549*theta_n.^5 - 37147.8947*theta_n.^6 + 38379.18127*theta_n.^7-16515.05308*theta_n.^8); % There is a typo in this eqn. in the original paper
	dudt_n = dudt_n./(1-48.09287227*theta_n+1017.234804*theta_n.^2-10481.80419*theta_n.^3+59431.3*theta_n.^4-195881.6488*theta_n.^5+...
		374577.3152*theta_n.^6 - 385821.1607*theta_n.^7 + 165705.8597*theta_n.^8);

	% Define the OCV for the negative electrode
	U_n   = 0.7222 + 0.1387*theta_n + 0.029*theta_n.^0.5 - 0.0172./theta_n + 0.0019./theta_n.^1.5 + 0.2808*exp(0.9-15*theta_n)-0.7984*exp(0.4465*theta_n - 0.4108);
	U_n   = U_n + (T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn)-param.Tref).*dudt_n;
end