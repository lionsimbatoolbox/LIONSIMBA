%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function is used to get a measurement of the SOC according to the
% internal states. This function assumes that all the states are
% measurable.

function Sout = internalSOCestimate(cs_average_t,param,i)
% Check if Fick's law of diffusion is used. This is required to define the
% correct way how to evaluate the SOC.
if(param{i}.SolidPhaseDiffusion~=3)
    cs_average = cs_average_t{i}(end,param{i}.Np+1:end);
else
    start_index = param{i}.Nr_p*param{i}.Np+1;
    end_index   = start_index+param{i}.Nr_n-1;
    cs_average  = zeros(param{i}.Nn,1);
    for n=1:param{i}.Nn
        cs_average(n)   = 1/param{i}.Rp_n*(param{i}.Rp_n/param{i}.Nr_n)*sum(cs_average_t{i}(end,start_index:end_index));
        start_index     = end_index + 1;
        end_index       = end_index + param{i}.Nr_n;
    end
end
Csout  = sum(cs_average);
Sout   = 100*(1/param{i}.len_n*(param{i}.len_n/(param{i}.Nn))*Csout/param{i}.cs_maxn);
end