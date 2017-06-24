%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rootFinder is used to determine the presence of an event during the
% resolution of the P2D model. This function is invoked at every
% integration instant, if the AppliedCurrent = 2 flag is set. In this
% particular case it has been used to determine discontinuities in the
% applied current profile. However, it can be used also for other purposes.
% Please refer to the IDA guide for a better understanding of how this
% script works.

function [g, flag, new_data] = rootFinder(t,y,yp,data)
% Initialize the redurned values. Please refer to the documentation of
% IDARootFn in order to investigate the structure of this function

g           = 1;
new_data    = [];
flag        = 0;


if(data.param{1}.AppliedCurrent == 2)
    % Evaluate the value of the applied current at time t
    t_value     = data.param{1}.getCurr(t,data.t0,data.tf,y,data.param,data.param{1}.extraData);
    
    if(isfield(data,'prevI'))
        % Evaluate the value of the applied current at time t+t*eps
        t_value   = data.param{1}.getCurr(t*(1+1e-5),data.t0,data.tf,y,data.param,data.param{1}.extraData);
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