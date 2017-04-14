%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2017: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIRSTORDERDERIVATIVEMATRIX precomputes the matrices used to implement numerical differentiation using 9 points.
function [derivativeMatrix,r8fdx]=firstOrderDerivativeMatrix(xl,xu,n)
%  Compute the spatial increment
dx=(xu-xl)/(n-1);
r8fdx=1./(40320.*dx);
nm4=n-4;

mid_block = zeros(n-8,n);

first_row 	= [-109584 	+322560 	-564480 	+752640 	-705600 	+451584 	-188160 	+46080 	-5040];
second_row 	= [-5040 	-64224 		+141120 	-141120 	+117600 	-70560 		+28224 		-6720 	+720];
third_row 	= [+720 	-11520 		-38304 		+80640 		-50400 		+26880 		-10080 		+2304 	-240];
fourth_row 	= [-240 	+2880 		-20160 		-18144 		+50400 		-20160 		+6720 		-1440 	+144];

i_th_row 	= [+144 	-1536 		+8064 		-32256 		+0 			+32256 		-8064 		+1536 	-144];

n_min_3_row 	= [-144 	+1440 		-6720 		+20160 		-50400 		+18144 		+20160 		-2880 	+240];
n_min_2_row 	= [+240 	-2304 		+10080 		-26880 		+50400 		-80640 		+38304 		+11520 	-720];
n_min_1_row 	= [-720 	+6720 		-28224 		+70560 		-117600 	+141120 	-141120 	+64224 	+5040];
n_min_0_row 	= [+5040 	-46080 		+188160  	-451584 	+705600 	-752640 	+564480 	-322560 +109584];

first_block = [first_row;second_row;third_row;fourth_row];
first_block = [first_block zeros(4,n-9)];

row_index = 1;
for i=5:nm4
    mid_block(row_index,row_index:row_index+8)=i_th_row;
    row_index = row_index+1;
end

last_block = [n_min_3_row;n_min_2_row;n_min_1_row;n_min_0_row];
last_block=[zeros(4,n-9) last_block];

derivativeMatrix = [first_block;mid_block;last_block];

end