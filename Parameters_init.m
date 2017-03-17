%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS_INIT defines the parameters used in the simulation
% environment. Modify the values here in in order to parametrize the cell
% variables or change simulator parameters.

function param = Parameters_init



%% Constants
% Faraday Constant  [C/mol]
param.F         = 96487;
% Gas constant      [J / (mol K)]
param.R         = 8.314;

%% Sections thickness    [m]
% Aluminium current collector
param.len_al= 10e-6;
% Positive Electrode
param.len_p = 80e-6;
% Separator
param.len_s = 25e-6;
% Negative Electrode
param.len_n = 88e-6;
% Carbon current collector
param.len_co= 10e-6;

%% Thermal conductivities [ W / (m K) ]

% Aluminium current collector
param.Lambda_al = 237;
% Positive electrode
param.Lambda_p  = 2.1;
% Separator
param.Lambda_s  = 0.16;
% Negative Electrode
param.Lambda_n  = 1.7;
% Carbon current collector
param.Lambda_co = 401;

%% Electrolyte diffusion coefficients [m^2 / s]

% Positive side
param.Dp = 7.5e-10;
% Separator
param.Ds = 7.5e-10;
% Negative Side
param.Dn = 7.5e-10;

%% Density [kg / m^3 ]
% Aluminium current collector
param.rho_al = 2700;
% Positive electrode
param.rho_p  = 2500;
% Separator
param.rho_s  = 1100;
% Negative electrode
param.rho_n  = 2500;
% Carbon current collector
param.rho_co = 8940;

%% Environment temperature [K]
param.Tref   = 298.15;

%% Specific heat coefficients [ J / (kg K) ]

% Aluminium current collector
param.Cpal   = 897;
% Positive Electrode
param.Cpp    = 700;
% Separator
param.Cps    = 700;
% Negative Electrode
param.Cpn    = 700;
% Carbon current collector
param.Cpco   = 385;

%% Solid phase conductivity [S/m]
param.sig_al = 3.55e7;
param.sig_co = 5.96e7;


%% Porisity indices

% Positive electrode
param.eps_p     = 0.385;
% Separator
param.eps_s     = 0.724;
% Negative electrode
param.eps_n     = 0.485;


%% Volume fraction
param.eps_fi    = [0.025;0;0.0326];

%% Bruggeman coefficients

% Positive electrode
param.brugg_p   = 4;
% Separator
param.brugg_s   = 4;
% Negative electrode
param.brugg_n   = 4;

%% Solid diffusion coefficients [m^2 / s]

% Positive electrode
param.Dps       = 1e-14;
% Negative electrode
param.Dns       = 3.9e-14;

%% Particle surface area [m^2 / m^3]
% Positive electrode
a_p       = 885000;
% Separator
a_s       = 0;
%Negative electrode
a_n       = 723600;
% Do not remove. It is used in the code
param.a_i       = [a_p;a_s;a_n];

%% Transference number - Not available for positive/negative electrode
param.tplus     = 0.364;

%% Reaction rate constants [ m ^ 2.5 / (mol^0.5 s ) ]
% Positive electrode
param.k_p       = 2.334e-11;
% Separator
param.k_s       = 0;
% Negative electrode
param.k_n       = 5.031e-11;

%% Heat exchange coefficient [W/m^2 K]
param.hcell     = 1;

%% Maximum concentration of Li-ions in the solid phase [ mol/m^3 ]
% Positive electrode
param.cs_maxp   = 51554;
% Separator
param.cs_maxs   = 0;
% Negative electrode
param.cs_maxn   = 30555;

% Solid particle radius [m] - It's equal for the positive and negative electrode
param.Rp_p     = 2e-6;
param.Rp_n     = 2e-6;


% Solid phase conductivities
param.sig    = [100;... % Positive electrode
    0;...   % Separator
    100     % Negative electrode
    ];

% Effective solid phase conductivities
param.sig_eff = param.sig.*(1 - [param.eps_p;param.eps_s;param.eps_n] - param.eps_fi);

%% Activation Energy for Temperature Dependent Solid Phase Diffusion [ J / mol ]

% Positive electrode
param.EaDps  = 5000;

% Negative Electrode
param.EaDns  = 5000;

%% Activation Energy for Temperature Dependent Reaction Constant [ J / mol ]

% Positive electrode
param.Eakip  = 5000;

% Negative electrode
param.Eakin  = 5000;


%% Initial conditions

% Electrolyte Li-ions initial concentration [mol/m^3]

param.ce_init = 1000;

% Initial temperature of the cell [K]
param.T_init = 298.15;


%% Simulator parameters

% Select the model used for approximate the solid phase diffusion
% Admitted values are:
%                       1 - Parabolic approximation (two parameters model)
%
%                       2 - Higher-order polynomial (three parameters
%                       model)
%
%                       3 - Full order model
%

param.SolidPhaseDiffusion = 1;

% Enable or disable temperature dynamics

param.TemperatureEnabled = 1;

% Integration step [s]
param.integrationStep = 0.5;

% Cutoff voltage [V]
param.CutoffVoltage = 2.5;

% Cutover voltage [V]
param.CutoverVoltage = 4.3;

% Cutoff SOC [%]
param.CutoffSOC = 0.9;

% Cutover SOC [%]
param.CutoverSOC = 90;

% Number of volumes in the aluminium current collector
param.Nal   = 10;

% Number of volumes in the positive electrode
param.Np    = 10;%21

% Number of volumes in the separator
param.Ns    = 10;

% Number of volumes in the negative electrode
param.Nn    = 10;

% Number of volumes in the carbon current collector
param.Nco   = 10;

% If the full diffusion model is selected, these two parameters define the
% number of discretization points inside the solid particles.

% Number of points for the cathode
param.Nr_p = 10;

% Number of points for the anode
param.Nr_n = 10;

% Initial concentration of Li-ions in the solid phase [mol/m^3]

% Positive electrode
param.cs_p_init = 25455;

% Negative electrode
param.cs_n_init = 26127;

% Enable or disable the scope in the matlab command line
param.Scope = 1;

% Enable or disable the printing of header information
param.PrintHeaderInfo = 1;

% Input current profile.
% Admitted values are:
%                       1 - Constant input current between the initial time
%                           and final time.
%
%                       2 - Variable profile described in the
%                           getInputCurrent script as a function of t
%
%                       3 - Potentiostatic charge. In this case the current
%                           is considered a variable and the charge is
%                           carried out at constant potential.

param.AppliedCurrent = 1;

%% External functions

% This field can be used as an extra structure and it is passed to all the
% external scripts.
param.extraData = [];

% Define the name of the external function that has to be called to compute
% the value of applied current. This function is called during the
% integration process.
param.CurrentFunction = @getInputCurrent;

% Define the name of the external function used to compute physical and
% transport properties of the materials during simulation. Please refer to
% the existing functions to get insight for custom implementations.

% Electrolyte diffusion coefficients
param.ElectrolyteDiffusionFunction          = @electrolyteDiffusionCoefficients;
% Electrolyte conductivity coefficients
param.ElectrolyteConductivityFunction       = @electrolyteConductivity;
% Open circuit potential
param.OpenCircuitPotentialFunction          = @openCircuitPotential;
% Solid phase diffusion coefficient
param.SolidDiffusionCoefficientsFunction    = @solidPhaseDiffusionCoefficients;
% Reaction rates
param.ReactionRatesFunction                 = @reactionRates;

% If a function handle is provided, the function is called after each
% integration step. All the states relative to the current integration step
% are provided besides extraData structure and timing information. See
% socEstimator.m as example
param.SOC_estimation_function = @socEstimator;

%% Potentiostatic section
% This value (if set together with the AppliedCurrent flag to 3) is used to control the battery in a potentiostatic
% fashion.
param.V_reference = 4;

%% Tolerances

% Integrator tolerances
param.AbsTol = 1e-6;
param.RelTol = 1e-6;

%% Ageing parameters (TESTING PURPOSES, BETA VERSION)

param.EnableAgeing = 0;

% Initial SEI resistance value [Ohm m^2]
param.R_SEI     = 0.01;
%Molar weight                               [kg/mol]
%ATTENTION: In Development of First Principles Capacity Fade Model
%for Li-Ion Cells, Ramadass et al. the measurement unit of M_p is wrong as
%well as the number itself. Please refer to Review of models for predicting
%the cycling performance of lithium ion batterise, Santhanagopalan et al.
param.M_n               = 73e-3;
% Admittance                                [S/m]
param.k_n_aging         = 3.79e-7;
% Side reaction current density             [A/m^2]
param.i_0_jside         = 0.80e-10;
% Open circuit voltage for side reaction    [V]
param.Uref_s            = 0.4;
% 1C current for the particular chemistry [A/m^2]
param.I1C 				= 29.5;
% Weight used in the aging dynamics
param.w 				= 2;


% Set to 1 if the user wants to use tha Jacobian matrix during
% calculations.
param.UseJacobian       = 1;

% This value, if set, represents the Jacobian function used from the
% integrator. It has to be a class 'Function' of the CasADi package. If
% provided, with UseJacobian=1, it will be used for speed up the
% integration process. If not provided, with UseJacobian=1, the code will
% compute the Jacobian on its own.
param.JacobianFunction = [];


% Type of the DAE system returned by LIONSIMBA batteryModel.m script. This function is under development.
% Admitted values are:
%                       1 - The equations are returned in an analytical
%                       form, written as implicit DAEs, i.e.
%
%                       x_dot - f(x,z) | Time differential equations
%                       z - g(x,z)     | Algebraic equations
%
%                       2 - The equations are returned in an analytical
%                       form, where time differential equations are written
%                       in an explicit form, i.e.
%
%                       f(x,z)         | Time differential equations
%                       z - g(x,z)     | Algebraic equations
%
param.daeFormulation = 1;

end
