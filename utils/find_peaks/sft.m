function yt=sft(y)
% Slow Fourier transform function
% yt = sft(y)
% Input
%    y - Time series
% Output
%    yt - Discrete Fourier transform
N = length(y);
twopiN = -2*pi*sqrt(-1)/N;
for k=0:N-1
 temp = exp(twopiN*(0:N-1)*k);
 yt(k+1) = sum(y .* temp);
end
return;