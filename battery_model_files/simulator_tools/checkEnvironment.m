%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECKENVIRONMENT checks if the required tools are present in the host
% system.

function checkEnvironment(param,num_arg)

% Check for input arguments
if(num_arg==0)
    error('Type help startSimulation')
end

if(num_arg == 3 && param{1}.AppliedCurrent == 1)
    error('For constant input currents, please provide the current density as a parameter to startSimulation')
end

% Check if SUNDIALS is installed correctly
disp('Checking for SUNDIALS availability...');
result = exist('IDAFree');
if(result~=2)
    error('The SUNDIALS software does not seem to be present in the search path or installed on your system.\nPlease go to:\nhttps://computation.llnl.gov/casc/sundials/main.html\nand follow the instruction for installing SUNDIALS in MATLAB.');
else
    disp('SUNDIALS Package FOUND!');
end

%Check for CasADi
disp('Checking for CasADi availability...');

try
    import casadi.*
catch e
    error('The SUNDIALS software does not seem to be present in the search path or installed on your system.\nPlease go to:\nhttps://github.com/casadi/casadi/wiki/InstallationInstructions\nand follow the instruction for installing CasADi in MATLAB.');
end

disp('CasADi Package FOUND!');


% For each parameter structure, check if all the fields are set.
for i=1:length(param)
    %Check if the parameters have been set correctly
    [result, missing] = checkBatteryParameters(param{i});
    if(result==0)
        disp('Warning, there are missing parameters in the param array.')
        disp('Here below there is the list of such parameters:')
        for jj=1:length(missing)
            disp(missing{jj});
        end
        error('Please fix the problem and restart the script')
    end
end

end