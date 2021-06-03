function s=bsmooth2(a,w)
%  bsmooth(a,w) smooths vector a by a boxcar (rectangular window) of width w
%  with tapered edge smooth.  T. C. O'Haver, 2008.
v=ones(1,w)/w;
S=conv(a,v);
startpoint=round((length(v) + 1)/2);
endpoint=round(length(a)+startpoint-1);
s=S(startpoint:endpoint);L=length(a);
% Now smooth the edges with progressively smaller boxcars(1)=(a(1)+a(2))./2;
for k=2:startpoint-1,
    s(k)=mean(a(1:(2*k-1)));
    s(L-k+1)=mean(a(L-2*k+2:L));
ends(L)=(a(L)+a(L-1))./2;