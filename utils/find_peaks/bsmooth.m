function s=bsmooth(a,w)
%  Convolution-based boxcar smooth
%  bsmooth(a,w) smooths vector a by a boxcar (rectangular window) of width w
%  T. C. O'Haver, 1988.
v=ones(1,w);S=conv(a,v);
startpoint=round((length(v) + 1)/2);
endpoint=round(length(a)+startpoint-1);
s=S(startpoint:endpoint) ./ sum(v);