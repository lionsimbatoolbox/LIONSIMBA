function Sout = socEstimator(t,t0,tf,states,extraData,param)
%	socEstimator estimates SOC of the cell.
%	This function can be called after each integration step (within the time-stepping loop).
%	It can be used to implement a custom SOC estimation by changing the function handle in
%	Parameters_init suitably. It is mandatory to respect the sign of the function and to return
%	a scalar value. The Sout value will be concatenated and it will be available in the result array
%	at the end of the simulation.

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

if(param.SolidPhaseDiffusion~=3)
    cs_average = states.cs_average(end,param.Np+1:end);
else
    start_index = param.Nr_p*param.Np+1;
    end_index   = start_index+param.Nr_n-1;
    cs_average  = zeros(param.Nn,1);
    for n=1:param.Nn
        cs_average(n)   = 1/param.Rp_n*(param.Rp_n/param.Nr_n)*sum(states.cs_average(end,start_index:end_index));
        start_index     = end_index + 1;
        end_index       = end_index + param.Nr_n;
    end
end
% Estimates the SOC according to the current information of the states.
Csout= sum(cs_average);
Sout = 100*(1/param.len_n*(param.len_n/(param.Nn))*Csout/param.cs_maxn);

end
