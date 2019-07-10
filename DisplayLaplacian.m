function Z = DisplayLaplacian(n1,n2);

M = size(n1);
M = M(1);
P = M-1;
h = 2/M;

D2X = 2*diag(ones(P^2,1))-diag(ones(P^2-1,1),1)-diag(ones(P^2-1,1),-1);
for i = P+1:P:P^2,
  D2X(i,i-1) = 0;
  D2X(i-1,i) = 0;
end
D2Y = 2*diag(ones(P^2,1))-diag(ones(P^2-P,1),P) - diag(ones(P^2-P,1),-P);
lap = D2X+D2Y;
longbv1 = cumsum([0 h*n1(1:M,1)']);
longbv2 = cumsum([0 h*n2(1,1:M)]);
longbv3 = cumsum([longbv1(1,M+1) h*n2(M,1:M)]);
longbv4 = cumsum([longbv2(1,M+1) h*n1(1:M,M)']);
  
Z = zeros(M+1,M+1);
Z(1:M+1,1) = longbv1';
Z(1,1:M+1) = longbv2;
Z(M+1,1:M+1) = longbv3;
Z(1:M+1,M+1) = longbv4';

border1 = (1:P);
border2 = (0:P-1)*P+1;
border3 = (1:P)*P;
border4 = ((P-1)*P+1:P^2); 
bordervec = zeros(1,P^2);
bv1 = longbv1(2:M);
bv2 = longbv2(2:M);
bv3 = longbv3(2:M);
bv4 = longbv4(2:M);

bordervec(border1) = bordervec(border1)+bv1;
bordervec(border2) = bordervec(border2)+bv2;
bordervec(border3) = bordervec(border3)+bv3;
bordervec(border4) = bordervec(border4)+bv4;
bordervec = bordervec(:);

  DNX = zeros(P^2,M^2);
block = zeros(P,M);
for i = 1:P,
  block(i,i) = 1;
  block(i,i+1) = -1;
end
for i = 1:P,
  DNX((i-1)*P+1:i*P,(i-1)*M+1:i*M) = block;
  DNX((i-1)*P+1:i*P,(i)*M+1:(i+1)*M) = block;
end
dnx = .5*DNX * n1(:).*h;

DNY = zeros(P^2,M^2);
block = zeros(P,M);
for i = 1:P,
  block(i,i) = 1;
  block(i,i+1) = 1;
end
for i = 1:P,
  DNY((i-1)*P+1:i*P,(i-1)*M+1:i*M) = block;
  DNY((i-1)*P+1:i*P,(i)*M+1:(i+1)*M) = -block;
end
dny = .5*DNY * n2(:).*h;
%figure(1)
%subplot(211);plot(dnx);
%subplot(212);plot(dny);
insideZ = lap \ (dnx+dny+bordervec);

Z(2:M,2:M) = reshape(insideZ,P,P);
