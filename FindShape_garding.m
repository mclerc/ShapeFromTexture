function [sigma,kt,kb,tau] = FindShape_garding(x);

a = x(1);
b = x(2);
c = x(3);
d = x(4);

sigma = atan(c);
kt = cos(sigma).*(a./c-2);
kb = d./c./cos(sigma);
tau = b./c;
