function [param,surface_area_per_face_for_given_layers,mass_cell,Cp_avg] = compute_lumped_mass_and_Cp_avg_for_given_layer_fcn(no_of_layers,param)
% compute_lumped_mass_and_Cp_avg_for_given_layer_fcn computes the cell's mass & Cp_avg. Returns them by appending their values to param struct (to be used in lumped thermal model calcs)

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

% Ignoring edge effects (but, can they be ignored when no. of layers is small?)
surface_area_per_face_for_given_layers = param.overall_surface_area_for_given_layers/no_of_layers;

vol_fraction_pos = (1 - param.eps_p - param.eps_fi(1)); % volume fraction occupied by the positive porous electrode material
vol_fraction_sep = (1 - param.eps_s - param.eps_fi(2)); % volume fraction occupied by the separator polymer material
vol_fraction_neg = (1 - param.eps_n - param.eps_fi(3)); % volume fraction occupied by the negative porous electrode material

volume_al = surface_area_per_face_for_given_layers*ceil(0.5*no_of_layers)*param.len_al;
volume_pos = surface_area_per_face_for_given_layers*param.len_p*vol_fraction_pos*no_of_layers;
volume_sep = surface_area_per_face_for_given_layers*param.len_s*vol_fraction_sep*no_of_layers;
volume_neg = surface_area_per_face_for_given_layers*param.len_n*vol_fraction_neg*no_of_layers;
volume_cu = surface_area_per_face_for_given_layers*ceil(0.5*(no_of_layers+1))*param.len_cu;
volume_LiPF6 = surface_area_per_face_for_given_layers*(param.len_p*param.eps_p + param.len_s*param.eps_s + param.len_n*param.eps_n)*no_of_layers;
volume_pouch = param.pouch_length*param.pouch_width*param.len_pouch;
% volume of tab material is ignored. Is this important?

mass_al = param.rho_al*volume_al;
mass_pos = param.rho_p*volume_pos;
mass_sep = param.rho_s*volume_sep;
mass_neg = param.rho_n*volume_neg;
mass_cu = param.rho_cu*volume_cu;
mass_LiPF6 = param.rho_LiPF6*volume_LiPF6;
mass_pouch = 2*param.rho_pouch*volume_pouch; % factor of two, for each of the top & bottom pouch covers
% mass of tab material is ignored. Is this important?

mass_cell = mass_al + mass_pos + mass_sep + mass_neg + mass_cu + mass_LiPF6 + mass_pouch; % [kg]
param.mass_cell = mass_cell; % appending this to param struct and return to calling script

% Cp of pouch is currently not available and not accounted for in
% calculation of Cp_avg. Also, note that tab material is not accounted for in Cp_avg calcs.
% Are these important?
Cp_avg = (param.Cpal*mass_al + param.Cpp*mass_pos + param.Cps*mass_sep + param.Cpn*mass_neg + param.Cpcu*mass_cu + param.CpLiPF6*mass_LiPF6 + param.Cppouch*mass_pouch)/param.mass_cell;    % Specific heat capacity, J/(kg K)
param.Cp_avg = Cp_avg;  % appending this to param struct and return to calling script

end