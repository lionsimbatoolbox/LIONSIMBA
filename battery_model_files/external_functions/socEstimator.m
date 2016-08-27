%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOCESTIMATOR defines the script which can be called after each
% integration step. It can be used to implement custom SOC estimation
% processes. In Parameters_init.m it is possible to change the function
% handle. It is mandatory to respect the sign of the function and to return
% a scalar value. The Sout value will be concatenated and it will be given
% in the result array after the end of the simulation.

function Sout = socEstimator(t,t0,tf,states,extraData,param)

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
    Sout = 100*(1/param.len_n*(param.len_n/(param.Nn))*Csout/param.cs_max(3));

end