%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function rhs = getCasadiFormulation(I,startParameters)

try
    test_lionsimba_folder
catch
    error('It seems that you did not add to the Matlab path the battery_model_files directory and the folders therein. Please fix this problem and restart the simulation.')
end

% Version of LIONSIMBA
version       = '1.022';

if(isempty(startParameters))
    % Load battery's parameters if not provided by the user
    param{1} = Parameters_init;
else
    % Use user provided parameters
    param = startParameters;
end

% Check the input current value.
if(param{1}.AppliedCurrent==1)
    if(~isreal(I) || isnan(I) || isinf(I) || isempty(I))
        if(~isa(I,'casadi.MX'))
            error('The input current provided is a complex value or NaN. Please check the value and restart.')
        end
    end
else
    error('getCasadiFormulation accepts only galvanostatic currents')
end



% Check for environmental tool availability
checkEnvironment(param,nargin);

% If enabled, print the header information
if(param{1}.PrintHeaderInfo==1)
    headerInfo(version)
end

% If everything is ok, let's start to simulate.
rhs = mainCore(I,param);
end



function result = mainCore(I,param)
import casadi.*


% Get the total number of cells that have to be simulated.
n_cells         = length(param);



% Check if more cells are simulated when potentiostatic conditions are
% required.
if(param{1}.AppliedCurrent==3 && n_cells~=1)
    clc
    error('!!!ERROR!!! -- Potentiostatic simulations are only possible with single cell packs -- !!!ERROR!!!')
end

% Switch among the selected operating modes defining suitable getCurr
% function. If multiple cells are required, the variable current profile
% or its constant value are retreived from the first element of the
% parameters structures. This is valid because, when multiple cells are
% connected in series, they are crossed by the same amount of current.
switch(param{1}.AppliedCurrent)
    case 2
        param{1}.getCurr = param{1}.CurrentFunction;
    otherwise
        param{1}.getCurr = @(t,t0,tf,extra)I;
end

% The Matlab function fsolve is used to initially solve the algebraic
% equations by keeping constant the differential ones to their initial
% values at t0.
opt_fsolve              = optimset;
opt_fsolve.Display      = 'off';
opt_fsolve.FunValCheck  = 'on';

n_diff          = zeros(n_cells,1);
n_alg           = zeros(n_cells,1);
start_x_index   = 1;
start_xp_index  = 1;


Y0  = [];
YP0 = [];
% Perform several checks over the cells

xvarsym_tot     = {};
xpvarsym_tot    = {};
zvarsym_tot     = {};
paramsym_tot    = {};
for i=1:n_cells
    if(~isfield(param{i},'symbolic_parameters'))
        error('Whet getCasadiFormulation is called, the parameters structure must contain the field symbolic_parameters in which all the symbolical parameters have to be stored')
    end
    % When Fick's law of diffusion is used, at least 10 discretization
    % points are required. Raise an error if this condition is not met.
    if((param{i}.Nr_p<10 || param{i}.Nr_n<10) && param{i}.SolidPhaseDiffusion==3)
        error('The number of discrete points for the paricles must be at least 10 in both cathode and anode.')
    end
    
    param{i}.Nsum      = param{i}.Np + param{i}.Ns + param{i}.Nn;
    param{i}.Nsum_nos  = param{i}.Np + param{i}.Nn;
    
    % Define the discretization steps.
    param{i}.deltax_al     = 1 / param{i}.Nal;
    param{i}.deltax_p      = 1 / param{i}.Np;
    param{i}.deltax_s      = 1 / param{i}.Ns;
    param{i}.deltax_n      = 1 / param{i}.Nn;
    param{i}.deltax_co     = 1 / param{i}.Nco;
    
    % Compute the indices where the differential and algebraic variables
    % are stored.
    param{i} = computeVariablesIndices(param{i});
    
    % Preallocate the differentiation matrices used for the solid phase
    % potential.
    param{i} = solidPhaseDifferentiationMatrices(param{i});
    
    % Preallocate the differentiation matrices used for the solid phase
    % diffusion in case of Fick's law.
    param{i} = solidPhaseDiffusionDifferentiationMatrices(param{i});
    
    % Init the value of the injected current.
    param{i}.I = param{1}.getCurr(0,0,10,param{1}.extraData);
    
    %% Initial conditions
    
    % Initial conditions for the differential states
    [cs_average_init, ce_init, T_init, film_init, Q_init, n_diff(i)] =   differentialInitialConditions(param{i});
    
    % Initial guesses for the algebraic variables
    [x0_alg, n_alg(i)] =   algebraicInitialConditions(param{i});
    
    % Store the number of differential and algebraic variables for each cell.
    param{i}.ndiff = n_diff(i);
    param{i}.nalg  = n_alg(i);
    
    % Solve the algebraic equations to find a set of semi-consistent initial
    % conditions for the algebraic equations. This will help the DAE solver as
    % a warm startup.
    xvarsym         = MX.sym(['xvar_cell',num2str(i)],n_diff(i));
    xpvarsym        = MX.sym(['xpvar_cell',num2str(i)],n_diff(i));
    xpvarsym_tot    = {xpvarsym_tot{:}, xpvarsym};
    zvarsym         = MX.sym(['zvar_cell',num2str(i)],n_alg(i));
    algebraic_function = Function('alg_fun',{zvarsym,xvarsym,param{i}.symbolic_parameters},{algebraicStates(zvarsym,xvarsym(param{i}.ce_indices),xvarsym(param{i}.cs_average_indices),xvarsym(param{i}.Q_indices),xvarsym(param{i}.T_indices),xvarsym(param{i}.film_indices),param{i})});
    
    algebraic_root{i} = rootfinder('algebraic_variables_root_function','newton',algebraic_function);
    
    xvarsym_tot     = {xvarsym_tot{:}, xvarsym};
    zvarsym_tot     = {zvarsym_tot{:}, zvarsym};
    
    paramsym_tot    = {paramsym_tot{:}, param{i}.symbolic_parameters};
    % Build the initial values array for the integrator
    Yt0 = [ce_init;cs_average_init;T_init;film_init;Q_init];
    
    
    Y0  = [Y0;Yt0];
    YP0 = [YP0;zeros(size(Yt0))];
    
    
    % The x_index variable will be used in the battery model file
    % for indexing purposes
    param{i}.x_index    = (start_x_index:n_diff(i)+n_alg(i)+start_x_index-1);
    param{i}.xp_index   = (start_xp_index:n_diff(i)+start_xp_index-1);
    
    % Update the starting x_index value for the (possible) next cell
    start_x_index       = n_diff(i)+n_alg(i)+start_x_index;
    start_xp_index      = n_diff(i)+start_xp_index;
    
end

% Define the structure to be passed to the residual function
dati.param  = param;
dati.t0     = 0;
dati.tf     = 10;

[dx_tot, ~, ~] = batteryModel(0,vertcat([xvarsym_tot{:};zvarsym_tot{:}]),vertcat([xpvarsym_tot{:}]),dati);


result.x_vars       = vertcat([xvarsym_tot{:}]);
result.z_vars       = vertcat([zvarsym_tot{:}]);
result.xp_vars      = vertcat([xpvarsym_tot{:}]);
result.param_vars   = vertcat([paramsym_tot{:}]);
result.ode          = dx_tot(1:n_diff(1));
result.ae           = dx_tot(n_diff(1)+1:end);
result.x0           = Y0;
result.z0           = x0_alg;
result.rootFunction = algebraic_root;
result.param        = param;

end