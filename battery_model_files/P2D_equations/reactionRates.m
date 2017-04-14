%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTIONRATES  Evaluates the reaction rates of the cathode 
%                       and anode of the battery. 
%                       The measurement unit is [m^2.5 / (m^0.5 s)]
%
%   You can modify the script to meet your particular requirements.

function [k_pT, k_nT] = reactionRates(T,param)

if(param.TemperatureEnabled==1)
    k_pT     = param.k_p*exp(-param.Eakip/param.R*(1./T(param.Nal+1:param.Nal+param.Np)-1/param.Tref));
else
    k_pT     = param.k_p;
end


if(param.TemperatureEnabled==1)
    k_nT     = param.k_n*exp(-param.Eakin/param.R*(1./T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn)-1/param.Tref));
else
    k_nT     = param.k_n;
end

end