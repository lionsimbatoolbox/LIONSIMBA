%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECKINITIALSTATES  Checks if the initial state structure provided to the
% software is suitable.

function [Y0_existence,YP0_existence,Y0,YP0] = checkInitialStates(initialState)

if(~isempty(initialState))
    Y_field_existence        = isfield(initialState,'Y');
    YP_field_existence       = isfield(initialState,'YP');
    
    if(Y_field_existence==0 || YP_field_existence==0)
        clc
        error('When defining initial states structure, please remember to set the Y and YP field')
    end
    
    if(Y_field_existence==1 && YP_field_existence==1)
        if(~isempty(initialState.Y) && ~isempty(initialState.Y))
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