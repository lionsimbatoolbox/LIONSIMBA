%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BATTERYMODEL evaluates the set of DAEs to compute the Li-ion cell model.

function [dx_tot, flag, new_data] = batteryModel(t,x_tot,xp_tot,dati)
% These two flags are not used but required by IDAs solver.
flag            = 0;
new_data        = [];
% Retreive data from the UserData field of IDAs
param_tot       = dati.param;
t0              = dati.t0;
tf              = dati.tf;
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
    
    % Due to the presence of series-connected cells, theoretically they are
    % all crossed by the same current. Novertheless, each cell can be
    % driven with its own current. This is useful for ABMS development
    % where charge balancing is required.
    param.I     = param_tot{1}.getCurr(t,t0,tf,x_tot,param_tot,param_tot{i}.extraData);
    % Retreive differential variables
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
    % App current
    I_app       = x(param.Iapp_indices);
    
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
    % Check for potentiostatic charge
    if(param.AppliedCurrent==3)
        param.I = I_app;
    end
    
    
    % Check the number of cells in order to provide the correct BCs for the
    % thermal dynamics (if the thermal dynamics are active). When thermal
    % dynamics are on, the heat flux between two adjacent cells has to be
    % kept continuous.
    if(param.TemperatureEnabled == 1)
        if(n_cells==1)
            % If only one cell it is present, then the cell will exchange heat
            % with the surrounding environment.
            param.T_BC_sx = param.hcell*(param.Tref-T(1))/(param.deltax_al*param.len_al);
            param.T_BC_dx = -param.hcell*(T(end)-param.Tref)/(param.deltax_co*param.len_co);
        else
            % If more than one cell is used, than the exchange of heat will be
            % done accordingly with the neighbor cells. Only the external
            % cathode and anode will dissipate heat with the surrounding
            % environment, while the inner cells will keep the heat flux
            % continuous.
            if(i==1)
                param.T_BC_sx = param.hcell*(param.Tref-T(1))/(param.deltax_al*param.len_al);
                
                % Find the deltax of the last volume of the current cell
                deltax_sx = param.len_co*param.deltax_co;
                % Find the deltax of the first volume of next cell
                deltax_dx = param_tot{i+1}.len_al*param_tot{i+1}.deltax_al;
                % Find the overall distance between the two adjacent nodes
                % in the two different cells.
                deltax_tot = (deltax_sx/2+deltax_dx/2);
                % Find the Harmonic Mean of the heat transfer coefficients
                beta = (deltax_sx/2)/deltax_tot;
                Lambda_star = (param.Lambda_co*param_tot{i+1}.Lambda_al)/(beta*param_tot{i+1}.Lambda_al + (1-beta)*param.Lambda_co);
                % Set the BC on the right side of the current cell. Note
                % that x_next_cell(param_tot{i+1}.T_indices(1)) refers to the
                % T(1) of the next cell
                x_next_cell   = x_tot(param_tot{i+1}.x_index);
                param.T_BC_dx = Lambda_star*((x_next_cell(param_tot{i+1}.T_indices(1))-T(end))/deltax_tot)/(param.deltax_co*param.len_co);
            elseif(i>1 && i<n_cells)
                % Find the deltax of the last volume of the previous cell
                deltax_sx = param_tot{i-1}.len_co*param_tot{i-1}.deltax_co;
                % Find the deltax of the first volume of the current cell
                deltax_dx = param.len_al*param.deltax_al;
                % Find the overall distance between the two adjacent nodes
                % in the two different cells.
                deltax_tot = (deltax_sx/2+deltax_dx/2);
                % Find the Harmonic Mean of the heat transfer coefficients
                beta = (deltax_sx/2)/deltax_tot;
                Lambda_star = (param.Lambda_al*param_tot{i-1}.Lambda_co)/(beta*param.Lambda_al + (1-beta)*param_tot{i-1}.Lambda_co);
                % Set the BC on the left side of the current cell
                x_prev_cell   = x_tot(param_tot{i-1}.x_index);
                param.T_BC_sx = -Lambda_star*((T(1)-x_prev_cell(param_tot{i-1}.T_indices(end)))/deltax_tot)/(param.deltax_al*param.len_al);
                
                % Find the deltax of the last volume of the current cell
                deltax_sx = param.len_co*param.deltax_co;
                % Find the deltax of the first volume of the next cell
                deltax_dx = param_tot{i+1}.len_al*param_tot{i+1}.deltax_al;
                % Find the overall distance between the two adjacent nodes
                % in the two different cells.
                deltax_tot = (deltax_sx/2+deltax_dx/2);
                % Find the Harmonic Mean of the heat transfer coefficients
                beta = (deltax_sx/2)/deltax_tot;
                Lambda_star = (param.Lambda_co*param_tot{i+1}.Lambda_al)/(beta*param_tot{i+1}.Lambda_al + (1-beta)*param.Lambda_co);
                % Set the BC on the right side of the current cell
                x_next_cell   = x_tot(param_tot{i+1}.x_index);
                param.T_BC_dx = Lambda_star*((x_next_cell(param_tot{i+1}.T_indices(1))-T(end))/deltax_tot)/(param.deltax_co*param.len_co);
            else
                % Find the deltax of the last volume of the previous cell
                deltax_sx = param_tot{i-1}.len_co*param_tot{i-1}.deltax_co;
                % Find the deltax of the first volume of the current cell
                deltax_dx = param.len_al*param.deltax_al;
                % Find the overall distance between the two adjacent nodes
                % in the two different cells.
                deltax_tot = (deltax_sx/2+deltax_dx/2);
                % Find the Harmonic Mean of the heat transfer coefficients
                beta = (deltax_sx/2)/deltax_tot;
                Lambda_star = (param.Lambda_al*param_tot{i-1}.Lambda_co)/(beta*param.Lambda_al + (1-beta)*param_tot{i-1}.Lambda_co);
                % Set the BC on the left side of the current cell
                x_prev_cell   = x_tot(param_tot{i-1}.x_index);
                param.T_BC_sx = -Lambda_star*((T(1)-x_prev_cell(param_tot{i-1}.T_indices(end)))/deltax_tot)/(param.deltax_al*param.len_al);
                param.T_BC_dx = -param.hcell*(T(end)-param.Tref)/(param.deltax_co*param.len_co);
            end
        end
    end
    % Build the array of initial conditions for the algebraic equations
    x0_alg      = [jflux;Phis;Phie;js;I_app];
    
    % Get the residuals from the algebraic equations
    [dxalg,Up,Un,dudt_p,dudt_n,Keff,J_S] = algebraicStates(x0_alg,ce,cs_barrato,Q,T,film,param);
    
    % Evaluate the residual for the electrolyte concentration
    [resCe,rhsCe] = electrolyteDiffusion(t,ce,dCe,jflux + [zeros(param.Np,1);J_S],T,param);
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
    if(isa(param.I,'casadi.SX') && param.EnableAgeing==1)
        % Film thickness
        [resDfilm, rhsDfilm] = SEI_layer(dFilm, J_S, param);
        
        % Use the if_else CasADi sentence to determine a switching behavior of the film thickness dynamics according to the sign of the applied current density.
        resDfilm = if_else(param.I>=0,resDfilm,zeros(size(resDfilm,1),1));
        rhsDfilm = if_else(param.I>=0,rhsDfilm,zeros(size(rhsDfilm,1),1));
    elseif(param.EnableAgeing==1 && param.I >=0)
		% If the applied current is a numerical quantity, then perform the regular computations
        [resDfilm, rhsDfilm] = SEI_layer(dFilm, J_S, param);
    else
        resDfilm    = dFilm;
        rhsDfilm    = zeros(length(resDfilm),1);
    end
    
    % Check if the temperature dynamics are enabled.
    if(param.TemperatureEnabled==1)
        % Evaluate the residuals for the temperature
        [resdT, rhsT] = thermalModel(ce,Phie,Phis,Keff,jflux+ [zeros(param.Np,1);J_S],T,dT,Up,Un,dudt_p,dudt_n,param);
    else
        % Return constant value for the temperature
        resdT   = dT;
        rhsT    = zeros(length(dT),1);
    end
    
    if(param_tot{1}.daeFormulation==1)
        % Return the residuals array
        dx_tot = [dx_tot;resCe;resCs;resdT;resDfilm;resdQ;dxalg];
    elseif(param_tot{1}.daeFormulation==2)
        % Return the RHSs of the equations
        dx_tot = [dx_tot;rhsCe;rhsCs;rhsT;rhsDfilm;rhsQ;dxalg];
    end
end
end
