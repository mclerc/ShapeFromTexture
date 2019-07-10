function [xn,yn,zn] = FindNormal(theta,sigma,x,y);
%
% From Garding :
%   N = p cos(sigma) - t sin(sigma) 
% 
% coordinates of t = (t1,t2,t3)
% satisfy t.p = 0 and t2/t1 = tan(theta)
for i = 1:prod(size(theta)),
  if theta(i) > pi/2,
    theta(i) = theta(i) - pi;
    sigma(i) = -sigma(i);
  end %if theta(i) >pi/2;
end %for i
t1 = 1./sqrt(1+tan(theta).^2+(x+y.*tan(theta)).^2);
t2 = t1.*tan(theta);
t3 = -sign(x+y.*tan(theta)).*sqrt(1-t1.^2-t2.^2);
p1 = x./sqrt(1+x.^2+y.^2);
p2 = y./sqrt(1+x.^2+y.^2);
p3 = 1./sqrt(1+x.^2+y.^2);
xn = cos(sigma) .* p1 - sin(sigma).* t1;
yn = cos(sigma) .* p2 - sin(sigma).* t2;
zn = cos(sigma) .* p3 - sin(sigma).* t3;

max(max(abs(t1.*p1 + t2.*p2 + t3.*p3)))


