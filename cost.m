function [f,x] = cost(theta,phi,m);

c = cos(theta);
s = sin(theta);
% we exchanged the middle two lines because of the way we order m
Jx = [ c.^3      -2.*c.^2.*s     2.*c.*s.^2       c.*s.^2   ;...
       c.^2.*s   -2.*c.*s.^2     -2.*c.^2.*s      s.^3      ;...
       c.^2.*s   c.*(c.^2-s.^2)  -s.*(c.^2-s.^2)  -c.^2.*s  ;...
       c.*s.^2   s.*(c.^2-s.^2)  c.*(c.^2-s.^2)   -c.*s.^2  ];

Jy = [ c.^2.*s   c.*(c.^2-s.^2)  -s.*(c.^2-s.^2)  -c.^2.*s  ;...
       c.*s.^2   s.*(c.^2-s.^2)  c.*(c.^2-s.^2)   -c.*s.^2  ;...
       c.*s.^2   2.*c.^2.*s      -2.*c.*s.^2       c.^3     ;...  
       s.^3      2.*c.*s.^2      2.*c.^2.*s        c.^2.*s];

for dir = 1:length(phi),
  J(4*(dir-1)+(1:4),:) = cos(phi(dir)).*Jx +sin(phi(dir)).*Jy; 
end %for dir;
x = inv(J' * J) * J' * m;
f = norm(J * x - m).^2;

% (c) Copyright Maureen Clerc,  March 1, 1999
