function s=bsmooth(a,w)
%  Convolution-based boxcar smooth

%  T. C. O'Haver, 1988.
v=ones(1,w);
startpoint=round((length(v) + 1)/2);
endpoint=round(length(a)+startpoint-1);
s=S(startpoint:endpoint) ./ sum(v);