%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECKBATTERYPARAMETERS  checks if the parameter structure is valid or not.

function [result, missing] = checkBatteryParameters(param)
parameters_list {1} = 'F';
parameters_list {end+1} = 'R';
parameters_list {end+1} = 'len_al';
parameters_list {end+1} = 'len_p';
parameters_list {end+1} = 'len_s';
parameters_list {end+1} = 'len_n';
parameters_list {end+1} = 'len_co';
parameters_list {end+1} = 'Lambda_al';
parameters_list {end+1} = 'Lambda_p';
parameters_list {end+1} = 'Lambda_s';
parameters_list {end+1} = 'Lambda_n';
parameters_list {end+1} = 'Lambda_co';
parameters_list {end+1} = 'Dp';
parameters_list {end+1} = 'Ds';
parameters_list {end+1} = 'Dn';
parameters_list {end+1} = 'rho_al';
parameters_list {end+1} = 'rho_p';
parameters_list {end+1} = 'rho_s';
parameters_list {end+1} = 'rho_n';
parameters_list {end+1} = 'rho_co';
parameters_list {end+1} = 'Tref';
parameters_list {end+1} = 'Cpal';
parameters_list {end+1} = 'Cpp';
parameters_list {end+1} = 'Cps';
parameters_list {end+1} = 'Cpn';
parameters_list {end+1} = 'Cpco';
parameters_list {end+1} = 'sig_al';
parameters_list {end+1} = 'sig_co';
parameters_list {end+1} = 'eps_p';
parameters_list {end+1} = 'eps_s';
parameters_list {end+1} = 'eps_n';
parameters_list {end+1} = 'eps_i';
parameters_list {end+1} = 'eps_fi';
parameters_list {end+1} = 'brugg_p';
parameters_list {end+1} = 'brugg_s';
parameters_list {end+1} = 'brugg_n';
parameters_list {end+1} = 'Dps';
parameters_list {end+1} = 'Dns';
parameters_list {end+1} = 'a_i';
parameters_list {end+1} = 'tplus';
parameters_list {end+1} = 'k_p';
parameters_list {end+1} = 'k_s';
parameters_list {end+1} = 'k_n';
parameters_list {end+1} = 'hcell';
parameters_list {end+1} = 'cs_max';
parameters_list {end+1} = 'Rp_p';
parameters_list {end+1} = 'Rp_n';
parameters_list {end+1} = 'sig';
parameters_list {end+1} = 'sig_eff';
parameters_list {end+1} = 'EaDps';
parameters_list {end+1} = 'EaDns';
parameters_list {end+1} = 'Eakip';
parameters_list {end+1} = 'Eakin';
parameters_list {end+1} = 'ce_init';
parameters_list {end+1} = 'T_init';
parameters_list {end+1} = 'SolidPhaseDiffusion';
parameters_list {end+1} = 'TemperatureEnabled';
parameters_list {end+1} = 'CutoffVoltage';
parameters_list {end+1} = 'CutoverVoltage';
parameters_list {end+1} = 'CutoffSOC';
parameters_list {end+1} = 'CutoverSOC';
parameters_list {end+1} = 'Nal';
parameters_list {end+1} = 'Np';
parameters_list {end+1} = 'Ns';
parameters_list {end+1} = 'Nn';
parameters_list {end+1} = 'Nco';
parameters_list {end+1} = 'Nr_p';
parameters_list {end+1} = 'Nr_n';
parameters_list {end+1} = 'cs_p_init';
parameters_list {end+1} = 'cs_n_init';
parameters_list {end+1} = 'Scope';
parameters_list {end+1} = 'PrintHeaderInfo';
parameters_list {end+1} = 'AppliedCurrent';
parameters_list {end+1} = 'extraData';
parameters_list {end+1} = 'CurrentFunction';
parameters_list {end+1} = 'ElectrolyteDiffusionFunction';
parameters_list {end+1} = 'ElectrolyteConductivityFunction';
parameters_list {end+1} = 'OpenCircuitPotentialFunction';
parameters_list {end+1} = 'SolidDiffusionCoefficientsFunction';
parameters_list {end+1} = 'ReactionRatesFunction';
parameters_list {end+1} = 'SOC_estimation_function';
parameters_list {end+1} = 'V_reference';
parameters_list {end+1} = 'AbsTol';
parameters_list {end+1} = 'RelTol';
parameters_list {end+1} = 'EnableAgeing';
parameters_list {end+1} = 'R_SEI';
parameters_list {end+1} = 'M_n';
parameters_list {end+1} = 'k_n_aging';
parameters_list {end+1} = 'i_0_jside';
parameters_list {end+1} = 'Uref_s';
parameters_list {end+1} = 'I1C';
parameters_list {end+1} = 'w';


k = 1;

for i=1:length(parameters_list)
    esito = isfield(param,parameters_list{i});
    if(esito==0)
        missing{k} = parameters_list{i};
        k = k+1;
    end
end
if(k==1)
    result = 1;
    missing = '';
else
    result = 0;
end
end