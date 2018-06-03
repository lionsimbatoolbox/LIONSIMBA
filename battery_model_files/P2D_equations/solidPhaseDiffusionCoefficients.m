function [Dps_eff, Dns_eff] = solidPhaseDiffusionCoefficients(T,param)
% solidPhaseDiffusionCoefficients evaluates diffusion coefficients of the solid phase [m^2 /s].
% The user may modify the script to meet specific requirements.

%   This file is part of the LIONSIMBA Toolbox
%
%	Official web-site: 	http://sisdin.unipv.it/labsisdin/lionsimba.php
% 	Official GitHUB: 	https://github.com/lionsimbatoolbox/LIONSIMBA
%
%   LIONSIMBA: A Matlab framework based on a finite volume model suitable for Li-ion battery design, simulation, and control
%   Copyright (C) 2016-2018 :Marcello Torchio, Lalo Magni, Davide Raimondo,
%                            University of Pavia, 27100, Pavia, Italy
%                            Bhushan Gopaluni, Univ. of British Columbia, 
%                            Vancouver, BC V6T 1Z3, Canada
%                            Richard D. Braatz, 
%                            Massachusetts Institute of Technology, 
%                            Cambridge, Massachusetts 02142, USA
%   
%   Main code contributors to LIONSIMBA 2.0:
%                           Ian Campbell, Krishnakumar Gopalakrishnan,
%                           Imperial college London, London, UK
%
%   LIONSIMBA is a free Matlab-based software distributed with an MIT
%   license.

if(param.TemperatureEnabled>=1)
    Dps_eff     = param.Dps*exp(-param.EaDps/param.R*(1./T(param.Nal+1:param.Nal+param.Np)-1/param.Tref));
else
    Dps_eff     = param.Dps*ones(param.Np,1);
end

if(param.TemperatureEnabled>=1)
    Dns_eff     = param.Dns*exp(-param.EaDns/param.R*(1./T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn)-1/param.Tref));
else
    Dns_eff     = param.Dns*ones(param.Nn,1);
end

end
