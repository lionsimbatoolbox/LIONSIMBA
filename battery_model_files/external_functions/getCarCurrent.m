%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIONSAMBA example script
% getCarCurrent : this scripts represents an example how piecewise inputs
% have to be defined to run properly with LIONSIMBA. This script returns the value of the input current as a function of
% the time, states and parameters

function I = getCarCurrent(t,t0,tf,y,param,extra)
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
