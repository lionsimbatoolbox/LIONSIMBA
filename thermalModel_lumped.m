function [res_dT, rhsT] = thermalModel_lumped(~,cs_star_avg,Phis,~,~,T,dT,Up,Un,dUdT_p,dUdT_n,param,I_density)
%   thermalModel_lumped evaluates the set of equations for (lumped) thermal dynamics.

%   This file  is part of the LIONSIMBA Toolbox
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

[Up_avg,dUdT_p_avg,Un_avg,dUdT_n_avg] = param.OpenCircuitPotentialFunction(cs_star_avg,T,param,1); %  the final argument is a dummy argument (currently under debug and may be removed in the future)

% Determines whether use linear interpolation or the value at the CV
% center.
switch param.edge_values
    case 1
        Up_pos_cc = Up_avg(1);
        Un_neg_cc = Un_avg(end);
    case 2
        Up_pos_cc = 1.5*Up_avg(1) - 0.5*Up_avg(2);
        Un_neg_cc = 1.5*Un_avg(end) - 0.5*Un_avg(end-1);
end

cell_avg_OCP = Up_pos_cc - Un_neg_cc;

% Determines whether use linear interpolation or the value at the CV
% center.
switch param.edge_values
    case 1
        Phis_pos_cc = Phis(1);
        Phis_neg_cc = Phis(end);
    case 2
        Phis_pos_cc = 1.5*Phis(1) - 0.5*Phis(2);
        Phis_neg_cc = 1.5*Phis(end) - 0.5*Phis(end-1);
end

V_cell = Phis_pos_cc - Phis_neg_cc;

cell_current = I_density*param.overall_surface_area_for_given_layers;
Q_polarisation = abs(cell_current)*abs(cell_avg_OCP - V_cell);	% shall always be positive

if param.lumped_thermal_version == 1 % No entropic heat generation
    mCpdTdt = -(param.h_lumped*param.tab_area)*(T - param.Tref)...
              + Q_polarisation;

elseif param.lumped_thermal_version == 2 % With entropic heat generation
    dUdT_p_avg = dUdT_p'*ones(length(dUdT_p),1)/param.Np;
    dUdT_n_avg = dUdT_n'*ones(length(dUdT_n),1)/param.Nn;

    sign_input_density = evaluate_sign_input_density(param); % evaluates sign of input current/power density as per operating mode

    if sign_input_density==1
        Q_rev = cell_current*T*(dUdT_p_avg-dUdT_n_avg); % Charging Entropic heat term
    else
        Q_rev = -cell_current*T*(dUdT_p_avg-dUdT_n_avg);  % Discharging Entropic heat term
    end

    mCpdTdt = -(param.h_lumped*param.tab_area)*(T - param.Tref)...
              + Q_polarisation + Q_rev;
else
    error('Incorrect value for choice of thermal model. Exiting..\n')
end

rhsT    = mCpdTdt/(param.mass_cell*param.Cp_avg);
res_dT  = (dT - rhsT);

end
