function exit_reason = checkSimulationStopConditions(n_cells, Phis_t, cs_bar_t, T_t, param)
% checkSimulationStopConditions checks if a particular stop condition for the simulation is met

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

exit_reason = 0;
% Check stop conditions for each cell
for i=1:n_cells
    Phis_pos_cc_t = 1.5*Phis_t{i}(end,1) - 0.5*Phis_t{i}(end,2);
    Phis_neg_cc_t = 1.5*Phis_t{i}(end,end) - 0.5*Phis_t{i}(end,end-1);
    voltage = Phis_pos_cc_t - Phis_neg_cc_t;
    Sout    = internalSOCestimate(cs_bar_t,param,i);
    max_layer_temperature = max(T_t{1}(end,:));
    % Break conditions.
    if(voltage<param{i}.CutoffVoltage)
        if param{i}.suppress_status_prints == 0
            fprintf('\nCell #%d  below its Cutoff voltage. Stopping ...\n',i);
        end
        if exit_reason==0
            exit_reason=[];
        end
        exit_reason = [exit_reason;1];
    end
    
    if(voltage>param{i}.CutoverVoltage)
        if param{i}.suppress_status_prints == 0
            fprintf('\nCell #%d  above its Cutover voltage. Stopping ...\n',i);
        end
        if exit_reason==0
            exit_reason=[];
        end
        exit_reason = [exit_reason;2];
    end
    
    if(Sout<param{i}.CutoffSOC)
        if param{i}.suppress_status_prints == 0
            fprintf('\nCell #%d  below its Cutoff SOC. Stopping ...\n',i);
        end
        if exit_reason==0
            exit_reason=[];
        end
        exit_reason = [exit_reason;3];
    end
    
    if(Sout>param{i}.CutoverSOC)
        if param{i}.suppress_status_prints == 0
            fprintf('\nCell #%d  above its Cutover SOC. Stopping ...\n',i);
        end
        if exit_reason==0
            exit_reason=[];
        end
        exit_reason = [exit_reason;4];
    end
    
    if(max_layer_temperature > 0.99*param{i}.Tmax)
        if param{i}.suppress_status_prints == 0
            fprintf('\nCell #%d  above its Maximum Permitted Temperature. Stopping ...\n',i);
        end
        if exit_reason==0
            exit_reason=[];
        end
        exit_reason = [exit_reason;5];
    end
    
    if(param{i}.enable_csneg_Saturation_limit == 1)
        if(param{i}.SolidPhaseDiffusion~=3)
            max_cs_surface_neg = max(cs_bar_t{i}(end,param{i}.Np+1:end));
        else
            max_cs_surface_neg = max(cs_bar_t{i}(end, (param{i}.Np*param{i}.Nr_p) +1:end));
        end
        
        if (max_cs_surface_neg > param{i}.cs_sat_thresh*param{i}.cs_neg_saturation)
            if param{i}.suppress_status_prints == 0
                fprintf('\n\nCell #%d  above its safe surface concentration (saturation threshold) Stopping ...\n\n',i);
            end
            if exit_reason==0
                exit_reason=[];
            end
            exit_reason = [exit_reason;6];
        end
    end
end
end
