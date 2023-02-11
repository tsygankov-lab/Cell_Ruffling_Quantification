function [coeff, yfft, yrec] = fourier_shape_SH2(ybnd,M)

nbnd = size(ybnd,1);
coeff = [];
yfft = [];

A = zeros(nbnd,2*M+1);
A(:,1) = 1;
%A(:,1) = 0;
for m=2:M+1
    A(:,m)   = cos((m-1)*(1:nbnd)/nbnd);
    A(:,m+M) = sin((m-1)*(1:nbnd)/nbnd);
end

yfft = A\ybnd;

coeff = sqrt(yfft(2:M+1,:).^2+yfft(M+2:end,:).^2);
coeff = coeff./sum(coeff,1);

yrec = A*yfft;

end