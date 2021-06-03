% Script that demonstrates the triangulation method of peak measurement
% Requires gaussian.m, findpeaksTplot.m, and findpeaksG.m in the path.
format short g
format compact
x=5:.05:15;
h1=.2+rand(1); % Peaks heights 
h2=.2+rand(1);
h3=.2+rand(1);
y1=h1.*gaussian(x,8,1); % First peak
y2=h2.*gaussian(x,10,1); % Second peak
y3=h3.*gaussian(x,13,1); % Third peak
y=y1+y2+y3+.01.*randn(size(x)); % Total signal plus white noise
clf;
plot(x,y,'.b')
SlopeThreshold=0.0003;
AmpThreshold=0.1;
smoothwidth=15;
peakgroup=15;
smoothtype=3;
plots=1;
Pt=findpeaksTplot(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype);
P=findpeaksG(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype);
title('Triangulation method using findpeaksT.m')
disp('True peak parameters')
TruePeakAreas=[trapz(x,y1) trapz(x,y2) trapz(x,y3)]';
disp('           Peak     Position     Height       Width       Area');
True=[1 8 h1 1 trapz(x,y1); 2 10 h2 1 trapz(x,y2); 3 13 h3 1 trapz(x,y3)];
disp(True)
disp('Results using using triangulation (findpeaksT)')
disp('           Peak     Position     Height       Width       Area');
disp(Pt)
disp('Results using using top-half Gaussian fitting (findpeaksG)')
disp('           Peak     Position     Height       Width       Area');
disp(P)
disp(' ')
