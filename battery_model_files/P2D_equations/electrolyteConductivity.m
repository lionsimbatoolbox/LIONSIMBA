%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELECTROLYTECONDUCTIVITY  Evaluates the conductivity coefficients for the
% electrolyte phase. The measurement unit is [S/m]
%
%   Keff = ELECTROLYTECONDUCTIVITY(ce,T,param) evaluates
%   the conductivity coefficients for the anode, separator and cathode electrolyte phase of the
%   battery. You can modify the script to meet your particular needs.
%
%   The conductivity coefficients can be evaluated in isothermal case
%   (param.TemperatureEnabled=0) or adiabatic case
%   (param.TemperatureEnabled=1).
%
%   You can modify the way that the conductivity coefficients are computed, as
%   function of electrolyte concentration and temperature. The main script
%   will pass also the param array.

function Keff = electrolyteConductivity(ce,T,param,batterySection)

if(param.TemperatureEnabled==1)
    switch(batterySection)
        case'p'
            Keff = param.eps_p^param.brugg_p *(1e-4*ce.*((-10.5+0.668*1e-3*ce+0.494*1e-6*ce.^2) +...
                (0.074  -1.78*1e-5*ce -8.86*1e-10*ce.^2).*T + (-6.96*1e-5+2.8*1e-8*ce).*T.^2).^2);
        case 's'
            Keff = param.eps_s^param.brugg_s *(1e-4*ce.*((-10.5+0.668*1e-3*ce+0.494*1e-6*ce.^2) +...
                (0.074  -1.78*1e-5*ce -8.86*1e-10*ce.^2).*T + (-6.96*1e-5+2.8*1e-8*ce).*T.^2).^2);
        case 'n'
            Keff = param.eps_n^param.brugg_n *(1e-4*ce.*((-10.5+0.668*1e-3*ce+0.494*1e-6*ce.^2) +...
                (0.074  -1.78*1e-5*ce -8.86*1e-10*ce.^2).*T + (-6.96*1e-5+2.8*1e-8*ce).*T.^2).^2);
    end
else
        switch(batterySection)
        case'p'
            Keff = param.eps_p^param.brugg_p *(4.1253*1e-2 + 5.007*1e-4*ce - 4.7212*1e-7*ce.^2 +1.5094*1e-10*ce.^3 -1.6018*1e-14*ce.^4);
        case 's'
            Keff = param.eps_s^param.brugg_s *(4.1253*1e-2 + 5.007*1e-4*ce - 4.7212*1e-7*ce.^2 +1.5094*1e-10*ce.^3 -1.6018*1e-14*ce.^4);
        case 'n'
            Keff = param.eps_n^param.brugg_n *(4.1253*1e-2 + 5.007*1e-4*ce - 4.7212*1e-7*ce.^2 +1.5094*1e-10*ce.^3 -1.6018*1e-14*ce.^4);
        end
end