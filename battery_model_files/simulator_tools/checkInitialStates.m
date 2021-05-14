function [Y0_existence,YP0_existence,Y0,YP0] = checkInitialStates(initialState)
% checkInitialStates checks if the initial state struct provided by the user is suitable.

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

if(~isempty(initialState))
    Y_field_existence        = isfield(initialState,'Y');
    YP_field_existence       = isfield(initialState,'YP');

    if(Y_field_existence==0 || YP_field_existence==0)
        clc
        error('When defining initial states structure, please remember to set the Y and YP field')
    end

    if(Y_field_existence==1 && YP_field_existence==1)
        if(~isempty(initialState.Y) && ~isempty(initialState.YP))
            Y0              = initialState.Y;
            YP0             = initialState.YP;
            Y0_existence    = 1;
            YP0_existence   = 1;
        else
            Y0              = [];
            YP0             = [];
            Y0_existence    = 0;
            YP0_existence   = 0;
        end
    end
else
    Y0_existence    = 0;
    YP0_existence   = 0;
    Y0              = [];
    YP0             = [];
end
end
