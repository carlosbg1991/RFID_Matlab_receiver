x=5:.1:65;
y=modelpeaks2(x,[1 5 5 5 5],[1 1 1 1 1],[10 20 30 40 50],[3 3 3 3 3],[0 -5 -10 -15 -20]);
clf;
plot(x,y,'.')
SlopeThreshold=0.001;
AmpThreshold=0.1;
smoothwidth=7;
peakgroup=31;
smoothtype=3;
plots=1;
Pt=findpeaksTplot(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype);
P=findpeaksG(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype);
title('Triangulation method using findpeaksT.m')
disp('Results using using findpeaksT')
disp('     Peak   Position    Height     Width     Area');
disp(Pt)
disp('Results using using findpeaksG')
disp(P)
disp(' ')
disp('% Accuracy of peak area')

100.*(P(1,5)- Pt(1,5))/P(1,5)



