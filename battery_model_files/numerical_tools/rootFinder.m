function [g, flag, new_data] = rootFinder(t,y,yp,data)
%	rootFinder is used to determine the presence of an event during the
%	resolution of the P2D model. This function is invoked at every
%	integration instant, if the OperatingMode = 2 flag is set. In this
%	particular case it has been used to determine discontinuities in the
%	applied current profile. However, it can be used also for other purposes.
%	Please refer to the IDA guide for a better understanding of how this
%	script works.

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

% Initialize the redurned values. Please refer to the documentation of
% IDARootFn in order to investigate the structure of this function

g           = 1;
new_data    = [];
flag        = 0;


if(data.param{1}.OperatingMode == 4)
    % Evaluate the value of the applied current at time t
    t_value     = data.param{1}.getCurrentDensity(t,data.t0,data.tf,y,data.param,data.param{1}.extraData);
    
    if(isfield(data,'prevI'))
        % Evaluate the value of the applied current at time t+t*eps
        t_value   = data.param{1}.getCurrentDensity(t*(1+1e-5),data.t0,data.tf,y,data.param,data.param{1}.extraData);
        if(norm(t-data.prevT)<1e-3 && norm(t_value-data.prevI)>5e-1)
            g = 0;
        elseif(t>data.prevT)
            g = 1;
        else
            g = -1;
        end
    else
        g = 1;
    end
    data.prevI  = t_value;
    data.prevT  = t;
    new_data    = data;
end
end