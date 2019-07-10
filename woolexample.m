N =256;
n = 2;% number of subimages to optimize the scales
      % for SolveDeltaA = 2^(2*n)
P = 2^n;
image= LoadPicture('pull_outside_3');
image = image(1:256,51:306);
imdft = fftshift(fft2(image.^2));
variance = abs(smooth(imdft,100));
perspim= image./sqrt(variance);

disp = 1;
if disp
 figure(1)
 imagesc(perspim); axis image; colormap gray
 xlabel('x','FontSize', 20); ylabel('y','FontSize', 20);
 handle_axis = gca;      % recupere le handle de l'axe actif
 set(handle_axis,'FontSize', 14)
end % if disp 

phi = [0,pi/2];
mX = -.01;
MX = .01;
mY = -.01;
MY = .01;
ndir = length(phi);
eps1 = (MX-mX)./N;
eps2 = (MY-mY)./N;
 
 theta = pi/6*[0 1 2 3 4 5];
scale2 = [ 2 3  1 1 ];
scale1 = [ 1 1  2 3 ];

DeltaA0 = SolveDeltaA_wool(perspim,n,phi,eps1,eps2,theta,scale1,scale2);
[x0,y0] = meshgrid(mX+eps1.*(1:N),mY+eps2.*(1:N));

Z_0 = MakeSurface(DeltaA0,mX,mY,eps1,eps2,phi,'');
 
beep

