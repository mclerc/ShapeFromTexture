function [result,mask] = gabor2d(imdft,s1,s2,theta,alpha,type,sigma2);
% IT IS NOT NORMALIZED (NORMALIZATION IS NOT NECESSARY IN OUR APPLICATION)
% CAUTION : THE ANGLE TURNS IN THE CLOCKWISE DIRECTION 
% (CONTRARY TO USUAL PRACTISE IN MATH)
% input : the fftshifted fft of the image

dims = size(imdft);
ctr = ceil((dims+0.5)/2);
switch type
  case 0, Name = 'normal';
  case 1, Name = 'derb11'; 
  case 2, Name = 'derb12';
  case 3, Name = 'derb21';
  case 4, Name = 'derb22';
  case 5, Name = 'derx';
  case 6, Name = 'dery';
end

xi = 2*pi;
ct = cos(theta);
st = sin(theta);
ca = cos(alpha);
sa = sin(alpha);
% coefficients of B = rot(-theta) * diag(s1,s2) * rot(-alpha)
b11 = s1.*ct.*ca - s2.*st.*sa;
b12 = s1.*ct.*sa + s2.*st.*ca;
b21 = -s1.*st.*ca - s2.*ct.*sa;
b22 = -s1.*st.*sa + s2.*ct.*ca;
[wx,wy] = meshgrid( ([1:dims(2)]-ctr(2))./(dims(2)/2), ...
    ([1:dims(1)]-ctr(1))./(dims(1)/2) );
% apply tranpose(B) to the vector (wx,wy)
rotwx = b11.*wx + b21.*wy;
rotwy = b12.*wx + b22.*wy;
analytic = rotwx >0;
mask = exp(-sigma2/2.*((rotwx-xi).^2+rotwy.^2));
%mask = mask.*analytic;

% calculation of the different derivatives
if strcmp(Name,'derb11'),
  mask = -sigma2.*wx.*(rotwx-xi).*mask;
elseif strcmp(Name,'derb12'),
  mask = -sigma2.*wx.*rotwy.*mask;
elseif strcmp(Name,'derb21'),
 mask = -sigma2.*wy.*(rotwx-xi).*mask;
elseif strcmp(Name,'derb22'),
 mask = -sigma2.*wy.*rotwy.*mask;
elseif strcmp(Name,'derx'),
 mask = -i.*pi.*wx.*mask;
elseif strcmp(Name,'dery'),
 mask = -i.*pi.*wy.*mask;
end

% inverse Fourier transform
    banddft = imdft.*conj(mask);
    result = ifft2(ifftshift(banddft));

% Copyright (c) Maureen Clerc, December 1998



