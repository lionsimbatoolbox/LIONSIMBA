% LIONSIMBA example script
% CC_CV_charge scenario: this script provides an example to run a
% constant current-constant voltage charge of a cell.

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

% Clear the workspace
clear
close all

%% Parameters
% Define the integration times.
t0 = 0;
tf = 10^6;

% Define the initial state structure
initialState.Y  = [];
initialState.YP = [];

% Define the parameters structure.
param{1}               = Parameters_init;

param{1}.Np            = 10;
param{1}.Ns            = 10;
param{1}.Nn            = 10;

param{1}.hcell         = 1;
param{1}.Tref          = 298.15;

param{1}.AbsTol        = 1e-6;
param{1}.RelTol        = 1e-6;

param{1}.CutoffSOC     = 20;

param{1}.SolidPhaseDiffusion = 1;

if(param{1}.SolidPhaseDiffusion == 3)
    multp = param{1}.Nr_p;
    multn = param{1}.Nr_n;
else
    multp = 1;
    multn = 1;
end

% Define 1C rate of current
I1C = 29.23;

% Run with a 1.5C-rate in the CC stage
C_rate = 1.5;
%% Discharge section FDM
% Discharge the battery to the 20% of SOC
out = startSimulation(t0,tf,initialState,-25,param);

% Store the Jacobian Matrix
param{1}.JacobianFunction = out.JacobianFun;

initialState = out.initialState;

% Redefine the cutoff SOC in order to avoid simulation interruptions
param{1}.CutoffSOC     = 2;

% Rest the battery with no current applied
out2 = startSimulation(0,5000,initialState,0,param);

initialState = out2.initialState;


delta_p1 = param{1}.cs_p_init - mean(out2.cs_surface{1}(end,1:param{1}.Np));
delta_n1 = param{1}.cs_n_init - mean(out2.cs_surface{1}(end,param{1}.Np+1:end));
ce_final1 = out2.ce{1}(end,end);

gigi{1} = Parameters_init(out2.SOC{1}(end));

error_p1 = gigi{1}.cs_p_init-mean(out2.cs_surface{1}(end,1:param{1}.Np));
error_n1 = gigi{1}.cs_n_init-mean(out2.cs_surface{1}(end,param{1}.Np+1:end));

param{1}.CutoverSOC = out.SOC{1}(1);

out3 = startSimulation(0,50000,initialState,25,param);
param{1}.CutoverSOC = out.SOC{1}(1)*1.1;
out4 = startSimulation(0,5000,out3.initialState,0,param);

delta_p2 = param{1}.cs_p_init - mean(out4.cs_surface{1}(end,1:gigi{1}.Np));
delta_n2 = param{1}.cs_n_init - mean(out4.cs_surface{1}(end,gigi{1}.Np+1:end));
ce_final2 = out4.ce{1}(end,end);


gigi{1} = Parameters_init(out4.SOC{1}(end));
error_p2 = gigi{1}.cs_p_init-mean(out4.cs_surface{1}(end,1:param{1}.Np));
error_n2 = gigi{1}.cs_n_init-mean(out4.cs_surface{1}(end,param{1}.Np+1:end));



