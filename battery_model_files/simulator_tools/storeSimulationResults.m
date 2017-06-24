%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% storeSimulationResults stores simulation results into internal structures.

function [ce_t,cs_barrato_t,T_t,jflux_t,Phis_t, Phie_t, cs_star_t, SOC_t, film_t, js_t,Up_t,Un_t,R_int_t,app_current_t,Voltage_t,SOC_estimated_t,Qrev_t,Qrxn_t,Qohm_t,Q_t,tot_voltage,dudtp_t, dudtn_t,t_tot] =...
    storeSimulationResults(n_cells,ce_t,cs_barrato_t,T_t,jflux_t,Phis_t, Phie_t, cs_star_t, SOC_t, film_t, js_t,app_current_t,Voltage_t,SOC_estimated_t,Up_t,Un_t,R_int_t,Qrev_t,Qrxn_t,Qohm_t, Q_t,dudtp_t,dudtn_t, t_tot, y, t,SOC_estimate,t0,tf, param)
tot_voltage = 0;
for i=1:n_cells
    [ce_t{i}, cs_barrato_t{i}, T_t{i}, jflux_t{i}, Phis_t{i}, Phie_t{i}, cs_star_t{i}, SOC_t{i}, film_t{i}, js_t{i},Q_t{i}, t_tot{i}] = retreiveData(ce_t{i}, cs_barrato_t{i}, T_t{i}, jflux_t{i},...
        Phis_t{i}, Phie_t{i}, cs_star_t{i}, SOC_t{i}, film_t{i}, js_t{i}, Q_t{i}, t_tot{i}, y(param{i}.x_index), t, param{i});
    % Estimate the SOC and evaluate the voltage using initial conditions
    % values.
    SOC_estimated_t{i}  = [SOC_estimated_t{i};SOC_estimate{i}(t,t0,tf,param{i},ce_t{i},cs_barrato_t{i},cs_star_t{i},Phie_t{i},Phis_t{i},jflux_t{i},T_t{i})];
    voltage             = (Phis_t{i}(end,1)-Phis_t{i}(end,end));
    tot_voltage         = tot_voltage + voltage;
    Voltage_t{i}        = [Voltage_t{i};voltage];
    app_current_t{i}    = [app_current_t{i};y(param{i}.x_index(end))];
    % Estimate the internal resistance
    [U_p,dudt_p,U_n,dudt_n]       = param{i}.OpenCircuitPotentialFunction(cs_star_t{i}(end,:)',T_t{i}(end,:)',param{i});
    Up_t{i}             = [Up_t{i};U_p'];
    dudtp_t{i}          = [dudtp_t{i};dudt_p'];
    dudtn_t{i}          = [dudtn_t{i};dudt_n'];
    Un_t{i}             = [Un_t{i};U_n'];
    % Get the OCV
    OCV                 = U_p(1)-U_n(end);
    R_int               = abs((voltage-OCV)/param{1}.getCurr(t,t0,tf,y,param,param{1}.extraData));
    R_int_t{i}          = [R_int_t{i};R_int];
    % Heat generation rates
    [Qrev, Qrxn, Qohm]  = heatGenerationRates(Phis_t{i}(end,:),Phie_t{i}(end,:),jflux_t{i}(end,:),T_t{i}(end,:),cs_star_t{i}(end,:),ce_t{i}(end,:),param{i});
    Qrev_t{i}           = [Qrev_t{i};Qrev];
    Qrxn_t{i}           = [Qrxn_t{i};Qrxn];
    Qohm_t{i}           = [Qohm_t{i};Qohm];
    
end
end

