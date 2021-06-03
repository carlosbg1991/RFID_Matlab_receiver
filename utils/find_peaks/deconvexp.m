function ydc=deconvexp(x,y,tc)
% function ydc=deconvgauss(x,y,tc) deconvolutes an exponential function of
% time constant 'tc' from vector y, returning the deconvoluted result
c=exp(-x./tc);    % exponential trailing convolution function, c
ydc=ifft(fft(y)./fft(c)).*sum(c);
