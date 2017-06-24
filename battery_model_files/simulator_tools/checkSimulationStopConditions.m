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
% This script is used to check several stop conditions for the battery

function exit_reason = checkSimulationStopConditions(n_cells, Phis_t, cs_bar_t, param )
exit_reason = 0;
% Check stop conditions for each cell
for i=1:n_cells
    voltage = Phis_t{i}(end,1)-Phis_t{i}(end,end);
    Sout    = internalSOCestimate(cs_bar_t,param,i);
    % Break conditions.
    if(voltage<param{i}.CutoffVoltage)
        disp(['Cell #',num2str(i),' below its Cutoff voltage. Stopping']);
        exit_reason = 1;
    end
    
    if(voltage>param{i}.CutoverVoltage)
        disp(['Cell #',num2str(i),' above its Cutover voltage. Stopping']);
        exit_reason = 2;
    end
    
    if(Sout<param{i}.CutoffSOC)
        disp(['Cell #',num2str(i),' below its Cutoff SOC. Stopping']);
        exit_reason = 3;
    end
    
    if(Sout>param{i}.CutoverSOC)
        disp(['Cell #',num2str(i),' above its Cutover SOC. Stopping']);
        exit_reason = 4;
    end
end
end