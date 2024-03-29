function y = rshift(x)
% rshift -- Circular right shift of 1-d signal
%  Usage
%    r = rshift(x)
%  Inputs
%    x   1-d signal
%  Outputs
%    r   1-d signal 
%        r(i) = x(i-1) except r(1) = x(n)
%

	n = length(x);
	y = [ x(n) x( 1: (n-1) )];

%
% Copyright (c) 1993. Iain M. Johnstone
%     
    
    
%   
% Part of WaveLab Version .701
% Built Tuesday, January 30, 1996 8:25:59 PM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@playfair.stanford.edu
%   
    
