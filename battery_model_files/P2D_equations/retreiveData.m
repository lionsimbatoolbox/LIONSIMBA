%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETREIVEDATA retreives the data which is returned in the results structure.

function [ce_t, cs_barrato_t, T_t, jflux_t, Phis_t, Phie_t, cs_star_t, SOC_t, film_t, js_t, Q_t,t_tot] = retreiveData(ce_t, cs_barrato_t, T_t, jflux_t, Phis_t, Phie_t,cs_star_t, SOC_t, film_t, js_t, Q_t,t_tot, y, t, param)
	% Extract differential variables after the integration process
	ce_t            = [ce_t;y(param.ce_indices)];
	cs_barrato_t    = [cs_barrato_t;y(param.cs_average_indices)];
	T_t             = [T_t;y(param.T_indices)];
	film_t          = [film_t;y(param.film_indices)];
    Q_t             = [Q_t;y(param.Q_indices)];
	% Extract the algebraic variables after the integration process
	jflux_t         = [jflux_t;y(param.jflux_indices)];
	Phis_t          = [Phis_t;y(param.Phis_indices)];
	Phie_t          = [Phie_t;y(param.Phie_indices)];
	js_t            = [js_t;y(param.js_indices)];

	% Check if Fick's law of diffusion is used. This is required to define the
    % correct way how to evaluate the SOC.
    if(param.SolidPhaseDiffusion~=3)
        cs_average = cs_barrato_t(end,param.Np+1:end);
    else
        start_index = param.Nr_p*param.Np+1;
        end_index   = start_index+param.Nr_n-1;
        cs_average  = zeros(param.Nn,1);
        for n=1:param.Nn
            cs_average(n)   = 1/param.Rp_n*(param.Rp_n/param.Nr_n)*sum(cs_barrato_t(end,start_index:end_index));
            start_index     = end_index + 1;
            end_index       = end_index + param.Nr_n;
        end
    end
    Csout           = sum(cs_average);
    Sout            = 100*(1/param.len_n*(param.len_n/(param.Nn))*Csout/param.cs_maxn);
    cs_star_t       = [cs_star_t;surfaceConcentration(cs_barrato_t(end,:)',jflux_t(end,:)',Q_t(end,:)',T_t(end,:)',param)'];
    SOC_t           = [SOC_t;Sout];
    t_tot           = [t_tot;t];
end