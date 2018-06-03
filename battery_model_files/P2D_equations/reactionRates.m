function [k_pT, k_nT] = reactionRates(T,param)
% REACTIONRATES  Reaction rates (k) of cathode and anode [m^2.5/(m^0.5 s)].
% The user may modify this script to meet specific requirements.

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
    k_pT     = param.k_p*exp(-param.Eakip/param.R*(1./T(param.Nal+1:param.Nal+param.Np)-1/param.Tref));
else
    k_pT     = param.k_p;
end

if(param.TemperatureEnabled>=1)
    k_nT     = param.k_n*exp(-param.Eakin/param.R*(1./T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn)-1/param.Tref));
else
    k_nT     = param.k_n;
end

end
