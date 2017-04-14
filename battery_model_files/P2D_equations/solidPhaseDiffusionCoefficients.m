%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLIDPHASEDIFFUSIONCOEFFICIENTS evaluates the diffusion coefficients for
% the solid phase. The measurement unit is [m^2 /s]
%
%   You can modify the script to meet your particular requirements.


function [Dps_eff, Dns_eff] = solidPhaseDiffusionCoefficients(T,param)

if(param.TemperatureEnabled==1)
    Dps_eff     = param.Dps*exp(-param.EaDps/param.R*(1./T(param.Nal+1:param.Nal+param.Np)-1/param.Tref));
else
    Dps_eff     = param.Dps*ones(param.Np,1);
end

if(param.TemperatureEnabled==1)
    Dns_eff     = param.Dns*exp(-param.EaDns/param.R*(1./T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn)-1/param.Tref));
else
    Dns_eff     = param.Dns*ones(param.Nn,1);
end

end