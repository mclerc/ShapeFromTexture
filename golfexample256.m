N = 256;
n = 3;% number of subimages to optimize the scales
      % for SolveDeltaA = 2^(2*n)
P = 2^n;
image= LoadPicture('golf');
image = image(1:2:2*N,1:2:2*N);
imdft = fftshift(fft2(image.^2));
% smooth : 1000 for the sea, 100 for sdbp1_128
variance = abs(smooth(imdft,400));
perspim= image./sqrt(variance);

figure(1)
imagesc(perspim); axis image; colormap gray
xlabel('x','FontSize', 20); ylabel('y','FontSize', 20);
handle_axis = gca;      % recupere le handle de l'axe actif
set(handle_axis,'FontSize', 14)

phi = [0, pi/2];
mX = -.01;
MX = .01;
mY = -.01;
MY = .01;
ndir = length(phi);
eps1 = (MX-mX)./N;
eps2 = (MY-mY)./N;
 
theta = [0 pi/6 pi/3 pi/2 2*pi/3 5*pi/6];
scale1 = [1.2 2   3.5 5];
scale2 = [1.2 1.5 2   3];

DeltaA0 = SolveDeltaA_golf(perspim,n,phi,eps1,eps2,theta,scale1,scale2);
[x0,y0] = meshgrid(mX+eps1.*(1:N),mY+eps2.*(1:N));

Z_0 = MakeSurface(DeltaA0,mX,mY,eps1,eps2,phi,'golf');
