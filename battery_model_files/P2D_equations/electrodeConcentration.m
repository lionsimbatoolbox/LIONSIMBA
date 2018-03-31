function [ddCs, rhsCs] = electrodeConcentration(dCs,cs_barrato,T,jflux,param)
% electrodeConcentration describes the ODE of the concentration of lithium ions within the
% electrodes.

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

if(param.SolidPhaseDiffusion==1 || param.SolidPhaseDiffusion==2)
    % Cathode
    rhsCs_p  =((-3/param.Rp_p)*jflux(1:param.Np));
    ddCs_p   = dCs(1:param.Np) - rhsCs_p;
    
    % Anode
    rhsCs_n  = ((-3/param.Rp_n)*jflux(param.Np+1:end));
    ddCs_n   = dCs(param.Np+1:end) - rhsCs_n;
else
    switch param.SolidPhaseDiffusionNumericalScheme
        % Use the FDM method for the solid phase diffusion
        case 1
            [rhsCs_p, rhsCs_n, ddCs_p, ddCs_n] = FDM9orderElectrodeDiffusion(T, cs_barrato, jflux, dCs, param);
        % Use the spectral method for the discretization of the solid phase
        % diffusion
        case 2
            [rhsCs_p, rhsCs_n, ddCs_p, ddCs_n] = spectralMethodElectrodeDiffusion(T, cs_barrato, jflux, dCs, param);
    end
end
rhsCs   = [rhsCs_p;rhsCs_n];
ddCs    = [ddCs_p;ddCs_n];
end
