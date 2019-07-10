function Z = MakeSurface(DeltaA,mX,mY,eps1,eps2,phi,type);

disp = 1;

if nargin<7,
	type = '';
end
[N,M,p,ndir] = size(DeltaA);
m = zeros(4*ndir,1);
zones = 32;
%  for izone = 1:zones, %   vertical
%ipix = (izone-1)*N/zones+1 : izone*N/zones;
%  for jzone = 1:zones, % horizontal
%jpix = (jzone-1)*N/zones+1 : jzone*N/zones;
  for izone = 1:zones-2, %   vertical
ipix = (izone)*N/zones+1 : (izone+1)*N/zones;
  for jzone = 1:zones-2, % horizontal
jpix = (jzone)*N/zones+1 : (jzone+1)*N/zones;
for dir = 1:ndir,
        m(4*(dir-1) + (1:4)) = mean(mean(DeltaA(ipix,jpix,:,dir)));
end % for dir
% to avoid falling into a local minimum
t1 = fminbnd('cost',0,pi-.0001,[],phi,m);
t2 = fminbnd('cost',pi/2,3*pi/2-.0001,[],phi,m);
      if cost(t1,phi,m) <= cost(t2,phi,m)
          theta(izone,jzone) = t1;
	else
	  theta(izone,jzone) = t2;
	end % if cost
	theta(izone,jzone) = rem(theta(izone,jzone),pi);
% get the 4 optimal parameters param 
% (corresponding to alpha, beta, gamma, delta in Gardings paper)
      [bogus,param] = cost(theta(izone,jzone),phi,m);       
% get the shape
      [sigma(izone,jzone),kt(izone,jzone),...
           kb(izone,jzone),tau(izone,jzone)] = FindShape_garding(param);
  end % for jzone
end % for izone
[x,y] = meshgrid(mX + eps1*N/zones.*((2:zones-1)-.5),...
				 mY + eps2*N/zones.*((2:zones-1)-.5));
[xn,yn,zn] = FindNormal(theta,sigma,x,y);
clear result;
result(:,:,1) = xn'./zn';
result(:,:,2) = yn'./zn';
if strcmp(type,'golf'),
    rad = max(max(x));
	n1 = -result(:,:,1).*(x.^2+y.^2 <1.5*(rad)^2);
	n2 = -result(:,:,2).*(x.^2+y.^2 <1.5*(rad)^2); 
elseif strcmp(type,'pcyl7'),
	n1 = -result(:,:,1).*(y-x<.012).*(y-x>-.012).*(y+x<.012).*(y+x>-.012);
	n2 = -result(:,:,2).*(y-x<.012).*(y-x>-.012).*(y+x<.012).*(y+x>-.012); 
else
	n1 = -result(:,:,1);
	n2 = -result(:,:,2); 
end
 P = zones-2;
 sigma2 = (P/16)^2;
 gauss = exp(-((1:P)-P/2).^2./sigma2);
 gauss = gauss./sum(gauss);
 gauss = fftshift(gauss);
 gauss = rshift(gauss);
 m1 = myconv(gauss, n1);
 m2 = myconv(gauss,n2);
if strcmp(type,'pcyl7'),
	m1 = m1.*(y-x<.012).*(y-x>-.012).*(y+x<.012).*(y+x>-.012);
	m2 = m2.*(y-x<.012).*(y-x>-.012).*(y+x<.012).*(y+x>-.012); 
end
if disp,
 figure(4)
 quiver(y,x,m1,fliplr(-m2)); axis image;
xlabel('x','FontSize',20); ylabel('y','FontSize',20);
end % if disp
Z = DisplayLaplacian(-m1,-m2);

if disp
 figure(6)
 meshc(flipud(Z'));
 xlabel('x','FontSize',20);ylabel('y','FontSize',20);zlabel('z','FontSize',20)
 handle_axis = gca;      % recupere le handle de l'axe actif
 set(handle_axis,'FontSize', 14)

end % if disp

if disp
  figure(3)
 imagesc(Z'); axis image; colormap gray
end % if disp

