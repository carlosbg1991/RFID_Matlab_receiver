function g = clippedlorentzian(x,position,width,m)
% clippedlorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix position and width scalar
% Value of function is clipped at (limited to a maximum value of)
%  m. Shape is Lorentzian when m=inf. T. C. O'Haver, 2014
% Example: clippedlorentzian([1 2 3],2,2,.9) gives result [0.5 .9 0.5]
g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);
for k=1:length(x),if g(k)>m,g(k)=m;end,end