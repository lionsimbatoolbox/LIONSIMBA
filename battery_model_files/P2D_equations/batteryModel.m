function [dx_tot, flag, new_data] = batteryModel(t,x_tot,xp_tot,ida_user_data)
% batteryModel returns the set of residuals of the complete battery model.

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

% These two flags are not used but required by IDA(s) solver.
flag            = 0;
new_data        = [];
% Retrieve data from the UserData field of IDA(s)
param_tot       = ida_user_data.param;
t0              = ida_user_data.t0;
tf              = ida_user_data.tf;
% Empty the total array of residuals.
dx_tot          = [];

% Find the number of cells in the stack
n_cells = length(param_tot);

% Evaluate over all the cells.
for i=1:n_cells
    % Associate to the "param" variable the parameters structure related to
    % the current cell.
    param       = param_tot{i};
    % From the overall number of variables, extract those related to the
    % current cell.
    x       = x_tot(param.x_index);
    xp      = xp_tot(param.xp_index);

    % All series-connected cells carry the same current.
    % Nevertheless, each cell can be driven with its own current. This
    % is useful for ABMS development wherein charge balancing is required.

    % Retrieve differential variables
    % Electrolyte concentration
    ce          = x(param.ce_indices);
    % Average solid phase concentration
    cs_barrato  = x(param.cs_average_indices);
    % Temperature
    T           = x(param.T_indices);
    % Film resistance
    film        = x(param.film_indices);
    % Averaged flux
    Q           = x(param.Q_indices);

    % Ionic flux
    jflux       = x(param.jflux_indices);
    % Solid phase potential
    Phis        = x(param.Phis_indices);
    % Electrolyte potential
    Phie        = x(param.Phie_indices);
    % Side reaction flux
    js          = x(param.js_indices);
    % current density
    I_density   = x(param.curr_dens_indices);

    % Electrolyte concentration derivative
    dCe         = xp(param.ce_indices,1);
    % Average surface concentration derivative
    dCs         = xp(param.cs_average_indices,1);
    % Temperature derivative
    dT          = xp(param.T_indices,1);
    % Film thickness
    dFilm       = xp(param.film_indices,1);
    % Average flux
    dQ          = xp(param.Q_indices,1);

    if param.OperatingMode==1 || param.OperatingMode==4
        param.I_density = param_tot{1}.getCurrentDensity(t,t0,tf,x_tot,param_tot,param_tot{i}.extraData); %param.I_density gets updated every time-step through this important line of code, i.e. fetch the value (either constant or from external function file) at time 't'
    elseif param.OperatingMode==2 || param.OperatingMode==5
        param.P_density = param_tot{1}.getPowerDensity(t,t0,tf,param_tot{i}.extraData);   %param.P_density gets updated every time-step through this important line of code, i.e. fetch the value (either constant or from external function file) at time 't'
    elseif param.OperatingMode==3    % Check if operating in potentiostatic charge mode
        param.I_density = I_density;
    else
        error('Not a valid operating mode');
    end

    % Check the number of cells in order to provide the correct BCs for the
    % thermal dynamics (if the thermal dynamics are active). When thermal
    % dynamics are on, the heat flux between two adjacent cells has to be kept continuous.
    if(param.TemperatureEnabled == 1)
        if(n_cells==1)
            % If only one cell it is present, then the cell will exchange heat
            % with the surrounding environment.
            param.T_BC_sx = param.hcell*(param.Tref-T(1))/(param.deltax_al*param.len_al);
            param.T_BC_dx = -param.hcell*(T(end)-param.Tref)/(param.deltax_cu*param.len_cu);
        else
            % If more than one cell is used, than the exchange of heat will be
            % done accordingly with the neighbor cells. Only the external
            % cathode and anode will dissipate heat with the surrounding
            % environment, while the inner cells will keep the heat flux
            % continuous.
            if(i==1)
                param.T_BC_sx = param.hcell*(param.Tref-T(1))/(param.deltax_al*param.len_al);

                % Find the deltax of the last volume of the current cell
                deltax_sx = param.len_cu*param.deltax_cu;
                % Find the deltax of the first volume of next cell
                deltax_dx = param_tot{i+1}.len_al*param_tot{i+1}.deltax_al;
                % Find the overall distance between the two adjacent nodes
                % in the two different cells.
                deltax_tot = (deltax_sx/2+deltax_dx/2);
                % Find the Harmonic Mean of the heat transfer coefficients
                beta = (deltax_sx/2)/deltax_tot;
                Lambda_star = (param.Lambda_cu*param_tot{i+1}.Lambda_al)/(beta*param_tot{i+1}.Lambda_al + (1-beta)*param.Lambda_cu);
                % Set the BC on the right side of the current cell. Note
                % that x_next_cell(param_tot{i+1}.T_indices(1)) refers to the
                % T(1) of the next cell
                x_next_cell   = x_tot(param_tot{i+1}.x_index);
                param.T_BC_dx = Lambda_star*((x_next_cell(param_tot{i+1}.T_indices(1))-T(end))/deltax_tot)/(param.deltax_cu*param.len_cu);
            elseif(i>1 && i<n_cells)
                % Find the deltax of the last volume of the previous cell
                deltax_sx = param_tot{i-1}.len_cu*param_tot{i-1}.deltax_cu;
                % Find the deltax of the first volume of the current cell
                deltax_dx = param.len_al*param.deltax_al;
                % Find the overall distance between the two adjacent nodes
                % in the two different cells.
                deltax_tot = (deltax_sx/2+deltax_dx/2);
                % Find the Harmonic Mean of the heat transfer coefficients
                beta = (deltax_sx/2)/deltax_tot;
                Lambda_star = (param.Lambda_al*param_tot{i-1}.Lambda_cu)/(beta*param.Lambda_al + (1-beta)*param_tot{i-1}.Lambda_cu);
                % Set the BC on the left side of the current cell
                x_prev_cell   = x_tot(param_tot{i-1}.x_index);
                param.T_BC_sx = -Lambda_star*((T(1)-x_prev_cell(param_tot{i-1}.T_indices(end)))/deltax_tot)/(param.deltax_al*param.len_al);

                % Find the deltax of the last volume of the current cell
                deltax_sx = param.len_cu*param.deltax_cu;
                % Find the deltax of the first volume of the next cell
                deltax_dx = param_tot{i+1}.len_al*param_tot{i+1}.deltax_al;
                % Find the overall distance between the two adjacent nodes
                % in the two different cells.
                deltax_tot = (deltax_sx/2+deltax_dx/2);
                % Find the Harmonic Mean of the heat transfer coefficients
                beta = (deltax_sx/2)/deltax_tot;
                Lambda_star = (param.Lambda_cu*param_tot{i+1}.Lambda_al)/(beta*param_tot{i+1}.Lambda_al + (1-beta)*param.Lambda_cu);
                % Set the BC on the right side of the current cell
                x_next_cell   = x_tot(param_tot{i+1}.x_index);
                param.T_BC_dx = Lambda_star*((x_next_cell(param_tot{i+1}.T_indices(1))-T(end))/deltax_tot)/(param.deltax_cu*param.len_cu);
            else
                % Find the deltax of the last volume of the previous cell
                deltax_sx = param_tot{i-1}.len_cu*param_tot{i-1}.deltax_cu;
                % Find the deltax of the first volume of the current cell
                deltax_dx = param.len_al*param.deltax_al;
                % Find the overall distance between the two adjacent nodes
                % in the two different cells.
                deltax_tot = (deltax_sx/2+deltax_dx/2);
                % Find the Harmonic Mean of the heat transfer coefficients
                beta = (deltax_sx/2)/deltax_tot;
                Lambda_star = (param.Lambda_al*param_tot{i-1}.Lambda_cu)/(beta*param.Lambda_al + (1-beta)*param_tot{i-1}.Lambda_cu);
                % Set the BC on the left side of the current cell
                x_prev_cell   = x_tot(param_tot{i-1}.x_index);
                param.T_BC_sx = -Lambda_star*((T(1)-x_prev_cell(param_tot{i-1}.T_indices(end)))/deltax_tot)/(param.deltax_al*param.len_al);
                param.T_BC_dx = -param.hcell*(T(end)-param.Tref)/(param.deltax_cu*param.len_cu);
            end
        end
    end

    % Build the array of initial conditions for the algebraic equations
    x0_alg      = [jflux;Phis;Phie;js;I_density];

    % Get the residuals from the algebraic equations
    [dxalg,Up,Un,dudt_p,dudt_n,Keff,J_S] = algebraicStates(x0_alg,ce,cs_barrato,Q,T,film,param);

    % Evaluate the residual for the electrolyte concentration
    [resCe,rhsCe] = electrolyteDiffusion(ce,dCe,jflux + [zeros(param.Np,1);J_S],T,param);
    % Evaluate the residual for the average surface concentration
    [resCs, rhsCs] = electrodeConcentration(dCs,cs_barrato,T,jflux,param);

    % Check the type of Solid Diffusion Approximation to be used.
    if(param.SolidPhaseDiffusion==2)
        [resdQ, rhsQ] = volumeAveragedConcentrationFlux(dQ,Q,jflux,T,param);
    else
        resdQ   = dQ;
        rhsQ    = zeros(length(dQ),1);
    end

    % Ageing effects take place only during charging processes (as assumed here).
    % It is necessary to switch between the case in which the applied current density is a numerical
    % quantity and in the case in which is a symbolical quantity.
    if(isa(I_density,'casadi.SX') && param.EnableAgeing==1)
        [resDfilm, rhsDfilm] = SEI_layer(dFilm, J_S, param);         % Film thickness

        % Use the if_else CasADi sentence to determine a switching behavior of the film thickness dynamics according to the sign of the applied current density.
        % resDfilm = if_else(I_density>=0,resDfilm,zeros(size(resDfilm,1),1));
        rhsDfilm = if_else(I_density>=0,rhsDfilm,zeros(size(rhsDfilm,1),1));
    elseif(param.EnableAgeing==1 && I_density >0) % If the applied current is a numerical quantity, then perform the regular computations

        [resDfilm, rhsDfilm] = SEI_layer(dFilm, J_S, param);   % Film thickness
    else
        resDfilm    = dFilm;
        rhsDfilm    = zeros(length(resDfilm),1);
    end

    % Check if the thermal dynamics are enabled.
    switch param.TemperatureEnabled
	case 1 % PDE-based thermal model
        % Evaluate the residuals for the temperature
        [resdT, rhsT] = thermalModel_pde(ce,Phie,Phis,Keff,jflux+ [zeros(param.Np,1);J_S],T,dT,Up,Un,dudt_p,dudt_n,x_tot(end),param);
	case 2 % Reduced order lumped thermal model
	    % Surface concentrations
        cs_star = surfaceConcentration(cs_barrato,jflux,Q,T,param); % column vector
        cs_star_avg_pos = cs_star(1:param.Np)'*ones(param.Np,1)/param.Np;        % (scalar) average surface conc in positive electrode
        cs_star_avg_neg = cs_star(param.Np+1:end)'*ones(param.Nn,1)/param.Nn;    % (scalar) average surface conc in negative electrode
        cs_star_avg = [cs_star_avg_pos*ones(param.Np,1);cs_star_avg_neg*ones(param.Nn,1)]; % replicate the average surface concentration suitably and concatenate

        % Evaluate the residuals vector for the temperature
        [resdT, rhsT] = thermalModel_lumped(ce,cs_star_avg,Phis,Keff,jflux + [zeros(param.Np,1);J_S],T,dT,Up,Un,dudt_p,dudt_n,param,I_density);
	otherwise % Thermal dynamics disabled
        % Return constant value for the temperature
        resdT   = dT;
        rhsT    = zeros(length(dT),1);
    end

    if(param_tot{1}.daeFormulation==1)
        % Return the residuals array
        dx_tot = [dx_tot;resCe;resCs;resdT;resDfilm;resdQ;dxalg];
    elseif(param_tot{1}.daeFormulation==2)
        % Return the RHS of the equations
        dx_tot = [dx_tot;rhsCe;rhsCs;rhsT;rhsDfilm;rhsQ;dxalg];
    end
end
end
