function resultingMatrix = blockDiagonalMatrix(param,varargin)
% blockDiagonalMatrix provides an interface for creating block diagonal matrices

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
  prev_class  = class(varargin{1});
  
  % Check if the matrices provided for building the block diagonal are all variables
  % of the same class
  for i=1:nargin-1
    if(isa(varargin{i},prev_class))
      prev_class=class(varargin{i});
    else
      error('All the inputs arguments must be of the same class')
    end
  end
  
  % When running the model using symbolic variables, the creation of the block diagonal has to be carried out
  % differently if Octave is used
  if(isa(varargin{1},'casadi.SX') || isa(varargin{1},'casadi.MX'))
    if(param.isMatlab)
      resultingMatrix = blkdiag(varargin{:});
    else
      size_tot                                                          = [0,0];  
      for i=1:nargin-1
        size_tot = size_tot + size(varargin{i});
      end
      resultingMatrix                                                             = SX.zeros(size_tot);
      start_position = [1,1];
      % Create the block diagonal
      for i=1:nargin-1
        resultingMatrix(start_position(1):start_position(1)+size(varargin{i},1)-1,start_position(2):start_position(2)+size(varargin{i},2)-1) = varargin{i};
        start_position(1) = start_position(1) + size(varargin{i},1);
        start_position(2) = start_position(2) + size(varargin{i},2);
      end
    end
  else
    resultingMatrix = blkdiag(varargin{:});
  end
    
end