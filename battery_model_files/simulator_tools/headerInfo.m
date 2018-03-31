function headerInfo(version)
% 	headerInfo prints the version information of LIONSIMBA

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

clc
fprintf('/------------------------------------\\\n')
fprintf('|                                    |\n')
fprintf('|            LIONSIMBA               |\n')
fprintf('|            Toolbox                 |\n')
fprintf('|            version %s           |\n',version)
fprintf('|                                    |\n')
fprintf('\\------------------------------------/\n')
fprintf('LIONSIMBA: A Matlab framework based on a finite volume model suitable for Li-ion battery design, simulation, and control\n')
fprintf('Copyright (C) 2015-%d by M.Torchio, L.Magni, B.Gopaluni, R.D.Braatz, and D.M.Raimondo\n',str2double(datestr(now,'yyyy'))+1)
fprintf('Send bug reports, questions or comments to davide.raimondo@unipv.it\n')
fprintf('Main code contributors to LIONSIMBA 2.0 Ian Campbell, Krishnakumar Gopalakrishnan, Imperial college London, London, UK')
fprintf('Updates available at the web page https://github.com/lionsimbatoolbox/LIONSIMBA\n\n')
end
