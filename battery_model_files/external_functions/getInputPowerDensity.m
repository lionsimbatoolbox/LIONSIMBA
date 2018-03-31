function P_density = getInputPowerDensity(t,t0,tf,x,param,extra)
%	getInputPowerDensity returns the value of the input current density as a function of time.
%
%       P = getInputPowerDensity(t,t0,tf,extra)
%
%       Inputs:
%               - t     : value of the current time step
%               - t0    : initial integration time
%               - tf    : final integration time
%               - extra : extra parameters
%       Outputs:
%               - P_density     : Applied power density [W/m^2]

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

% Define your own linear/nonlinear function of time for the applied power density

% P = (t-t0)/(tf-t0) *(-30) + 0;
% P_density = -90*sin(t+100/100);

% P_density = -10*sin(omega*t+phi);
frequency = 0.1;
P_density = -300*sin((2*pi*frequency)*t+(1/100));
if P_density == 0
    P_density = 1e-2;
end
    
end
