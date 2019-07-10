function [result,mask] = smooth(imdft,sigma2);
% IT IS NOT NORMALIZED (NORMALIZATION IS NOT NECESSARY IN OUR APPLICATION)
% CAUTION : THE ANGLE TURNS IN THE CLOCKWISE DIRECTION 
% (CONTRARY TO USUAL PRACTISE IN MATH)
% input : the fftshifted fft of the image

dims = size(imdft);
ctr = ceil((dims+0.5)/2);
[wx,wy] = meshgrid(([1:dims(2)]-ctr(2))./(dims(2)/2), ...
    ([1:dims(1)]-ctr(1))./(dims(1)/2) );
%analytic = wx >0;
mask = exp(-sigma2/2.*(wx.^2+wy.^2));
%mask = mask.*analytic;

% inverse Fourier transform
    banddft = imdft.*conj(mask);
    result = ifft2(ifftshift(banddft));

% Copyright (c) Maureen Clerc, December 1998



