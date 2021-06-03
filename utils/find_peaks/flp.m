function fy=flp(y,n)
% function fy=flp(y,n)
% Super-simple Fourier low pass filter, n=number of frequencies to pass.
% Returns fy, the simplification of y reconstructed from only the first n
% frequencies. This is intended mainly as an instructional tool, to
% demonstrate how to handle the format of Matlab's fft function.
%
% Example
%  x=0:.1:100;
%  y=sin(x)+sin(2*x);
%  subplot(2,1,1);plot(x,y);
%  subplot(2,1,2);plot(x,flp(y,20))
%
ffty=fft(y); % ffty is the fft of y
lfft=length(ffty); % Length of the FFT
ffty(n:lfft-n)=0; % All other frequencies set to zero, in both halves
fy=real(ifft(ffty)); % Real part of the inverse Fourier transform (ifft)