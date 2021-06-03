function g = clippedgaussian(x,pos,wid,m)
%  clippedgaussian(x,pos,wid.m) = clipped Gaussian centered on x=pos,
%  half-width=wid x may be scalar, vector, or matrix, pos and wid both
%  scalar. Value of function is clipped at (limited to a maximum value of)
%  m). Shape is Gaussian when m=inf. T. C. O'Haver, 2014
% Example: clippedgaussian([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
g = exp(-((x-pos)./(0.6005612.*wid)) .^2);
for k=1:length(x),if g(k)>m,g(k)=m;end,end