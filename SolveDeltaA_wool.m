function DeltaA= SolveDeltaA_cylge5(im,n,phi,eps1,eps2,theta,scale1,scale2);
% we multiply both directions by rho in this version
% we have simplified the partial derivatives with respect to
% s11,s12,s21,s22. They are calculated more directly than in 
% SolveDeltaA_gaborrho. 
% We still use rho in both directions.
% compared to SolveDeltaA_simple, we try to optimize the scale
% locally.
% Inputs:
%    im  : a N times N image
%     n  : number of subimages = 2^(2*n)
ndir = length(phi);
N = size(im);
N = N(1);
sigma2 = (N/3)^2;
gauss = exp(-((1:N)-N/2).^2./sigma2);
gauss = gauss./sum(gauss);
gauss = fftshift(gauss);
gauss = rshift(gauss);
imdft = fftshift(fft2(im));
%P = 6; % border safety
P = 20; % border safety 
%--------------- INITIALIZATION ---------------------
% scales to consider
%sc = [4 8 16 32 64];
sc = 32;
% works for :
%sc = 32;
nsc = length(sc);

sig2 = .5;
whichp = zeros(nsc,2^(2*n));
iscales = [];
a = 0;
% warplet transform over all scales
if 0,
for isc = 1:nsc,
  scale = sc(isc);
 s1 = scale*2;
 s2 = scale*2;
 for  it = 1:length(theta),
 WX0(:,:,isc,it)   = gabor2d(imdft,s1,s2,theta(it),a,0,sig2);
end;
end; % for isc

% local (in space) search for best scale
for i = 0 : 2^n-1, % vertical
  iint =i.*2.^(-n).*N+1 : (i+1) .* 2.^(-n).*N;
       for j = 0:2^n-1, % horizontal
jint = j.*2.^(-n).*N+1 : (j+1) .* 2.^(-n).*N;
	 p = 2^n*i+j+1;
       for isc = 1:nsc,
	 % keep only the best direction
	   WX = squeeze(abs(WX0(:,:,isc,:)).^2);
	     energy(isc) = max(mean(mean(WX(iint,jint,:))));
	end;
	iscale = keep(energy,1)
	iscales = union(iscales,iscale);	
	bogus = setdiff(union(whichp(iscale,:),p),[0]);
	whichp(iscale,1:length(bogus)) = bogus;
  end; % for j
end; % for i
end % if 0

iscales = [1];
whichp = (1:2^(2*n));
whichp

for iscale = iscales,
is = 0;
for it = 1:length(theta),
  t =theta(it);
  for isc = 1:length(scale1),
  s1 = sc(iscale).*scale1(isc);
  s2 = sc(iscale).*scale2(isc);
is = is + 1;

  b11 =  s1.*cos(t).*cos(a) - s2.*sin(t).*sin(a);
b12 =  s1.*cos(t).*sin(a) + s2.*sin(t).*cos(a);
b21 = -s1.*sin(t).*cos(a) - s2.*cos(t).*sin(a);
b22 = -s1.*sin(t).*sin(a) + s2.*cos(t).*cos(a);
WX0    = gabor2d(imdft,s1,s2,t,a,0,sig2);
WX11  = gabor2d(imdft,s1,s2,t,a,1,sig2);
WX12  = gabor2d(imdft,s1,s2,t,a,2,sig2);
WX21  = gabor2d(imdft,s1,s2,t,a,3,sig2);
WX22  = gabor2d(imdft,s1,s2,t,a,4,sig2);
dxWX  = gabor2d(imdft,s1,s2,t,a,5,sig2);
dyWX  = gabor2d(imdft,s1,s2,t,a,6,sig2);

WX11  = 2.*myconv(gauss,real(WX11.*conj(WX0)),P);
WX12  = 2.*myconv(gauss,real(WX12.*conj(WX0)),P);
WX21  = 2.*myconv(gauss,real(WX21.*conj(WX0)),P);
WX22  = 2.*myconv(gauss,real(WX22.*conj(WX0)),P);
dxWX  = 2.*myconv(gauss,real(dxWX.*conj(WX0)),P);
dyWX  = 2.*myconv(gauss,real(dyWX.*conj(WX0)),P);


for p = setdiff(whichp(iscale,:),[0]),
  ip = floor((p-1)./2^n);
jp = mod(p-1,2^n);
iint = ip.*2.^(-n).*N+1 : (ip+1) .* 2.^(-n).*N;
jint = jp.*2.^(-n).*N+1 : (jp+1) .* 2.^(-n).*N;

for dir = 1:ndir,
  vector(dir, p, is,:,:) = cos(phi(dir)).*dxWX(iint,jint)./eps1+sin(phi(dir)).*dyWX(iint,jint)./eps2;
end%for dir

matrix(p,is,1,:,:) = b11.*WX11(iint,jint) + b12.*WX12(iint,jint);
matrix(p,is,2,:,:) = b11.*WX21(iint,jint) + b12.*WX22(iint,jint);
matrix(p,is,3,:,:) = b21.*WX11(iint,jint) + b22.*WX12(iint,jint);
matrix(p,is,4,:,:) = b21.*WX21(iint,jint) + b22.*WX22(iint,jint);

end %for p
end % for isc
end % for it
end % for iscale


for  ipix = 1:N,
  for jpix = 1:N,
   ip = floor((ipix-1)./(2^(-n)*N));
   jp = floor((jpix-1)./(2^(-n)*N));
   p = 2^n*ip + jp + 1;
   i = ipix - ip.*2.^(-n).*N;
   j = jpix - jp.*2.^(-n).*N;
    for dir = 1:ndir,
	DeltaA(ipix,jpix,:,dir) = squeeze(matrix(p,:,:,i,j)) ...
		\squeeze(vector(dir,p,:,i,j));
    end % for dir
  end % for jpix
end % for ipix

% Copyright (c) Maureen Clerc, June 1999


