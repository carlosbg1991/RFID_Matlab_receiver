function ydc=deconvgauss(x,y,w)
% function ydc=deconvgauss(x,y,tc) deconvolutes a Gaussian function of
% width 'w' from vector y, returning the deconvoluted result
c=gaussian(x,0,w)+gaussian(x,max(x),w);   % Gaussian convolution function, c
ydc=ifft(fft(y)./fft(c)).*sum(c);

function g = gaussian(x,pos,wid)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.60056120439323.*wid)) .^2);