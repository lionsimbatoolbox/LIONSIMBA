%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function headerInfo(version)
clc
fprintf('/------------------------------------\\\n')
fprintf('|                                    |\n')
fprintf('|            LIONSIMBA               |\n')
fprintf('|            Toolbox                 |\n')
fprintf('|            version %s           |\n',version)
fprintf('|                                    |\n')
fprintf('\\------------------------------------/\n')
fprintf('Copyright (C) 2015-%d by M.Torchio, L.Magni, B.Gopaluni, R.D.Braatz and D.M.Raimondo\n\n',str2double(datestr(now,'yyyy'))+1)
fprintf('Send bug reports, questions or comments to marcello.torchio01@ateneopv.it\n')
fprintf('Updates available at the web page http://sisdin.unipv.it/labsisdin/lionsimba.php\n\n')
end