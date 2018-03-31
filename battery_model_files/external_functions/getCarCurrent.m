function I = getCarCurrent(t,t0,tf,x,param,extra)
%	getCarCurrent : this scripts represents an example how piecewise inputs
%	have to be defined to run properly with LIONSIMBA. This script returns the value of the input current as a function of
%	the time, states and parameters

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

if(t>=0 && t<=50)
    I = -29.5;
elseif(t>50 && t<=60)
    I = 14.75;
elseif(t>60 && t<=210)
    I = -14.75;
elseif(t>210 && t<=410)
    I = -29.5;
elseif(t>410 && t<=415)
    I = -58;
elseif(t>415 && t<=615)
    I = -37.5;
elseif(t>615)
    I = 14.75;
end
end