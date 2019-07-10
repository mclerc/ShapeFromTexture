function result = myconv(vector,image,P);
% separable convolution with FFT
if nargin < 3, 
  P = 0;
end
N = size(image);
N = N(1);
enveloppe = ones(1,N);
enveloppe(N-P+1:N) = (1+cos((1:P)./P*pi))./2;
enveloppe(1:P) = reverse(enveloppe(N-P+1:N));
enveloppe = enveloppe'*enveloppe;
smim = image.*enveloppe;
%smim = image;
N = length(vector);
% matrix should be N times N
filter = fft(vector);
filter = filter(:);
filtermat = filter * ones(1,N);
% first direction
bogus = ifft(fft(smim).*filtermat);
% second direction
result = ifft(fft(bogus').*filtermat)';
result = real(result);
%image(P+1:N-P,P+1:N-P) = result(P+1:N-P,P+1:N-P);
%result = image;

