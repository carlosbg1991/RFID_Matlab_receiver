function s=tsmooth(Y,w)
%  tsmooth(Y,w) smooths vector Y by a triangular function of halfwidth w
%  T. C. O'Haver, 1988.
v=ones(1,w);v=conv(v,v);
S=conv(Y,v);
startpoint=(length(v) + 1)/2;
endpoint=length(Y)+startpoint-1;
s=S(startpoint:endpoint) ./ sum(v);