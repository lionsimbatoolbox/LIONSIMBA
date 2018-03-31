function sign_input_density = evaluate_sign_input_density(param)
% evaluate_sign_input_density returns the sign of input current/power density as per operating mode

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

if param.OperatingMode==1 || param.OperatingMode==4
    sign_input_density = sign(param.I_density);
elseif param.OperatingMode==2 || param.OperatingMode==5
    sign_input_density = sign(param.P_density);
elseif param.OperatingMode==3
    sign_input_density = 1; % dummy value for CV mode (since usually we have CV charging)
else
    error('Not a valid operating mode.');
end

end
