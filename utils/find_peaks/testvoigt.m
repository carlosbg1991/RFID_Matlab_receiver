   % Demonstration of fitting two strongly overlapping Voigt peaks with
   % different alphas (ratio of Gaussian width to Lorentzian width), using
   % the unconstrained Voigt peak shape number 30. The peak positions are
   % recovered accurately but the alphas are not so precise because of
   % propagation of errors.
   format compact
   format short g
   x=0:.01:2;
   pos1=.8;
   pos2=1;
   WidthA=.3;
   WidthB=.4;
   ActuaAlpha1=WidthA/WidthB
   ActuaAlpha2=WidthB/WidthA
   hw=halfwidth(x,y); % Halfwidth of largest peak (to use as "extra")
   y=voigt(x,pos1,WidthA,WidthB)+.8.*voigt(x,pos2,WidthB,WidthA)+.01.*randn(size(x));
   % For a more stable fit in this case, use NumTrials=10 and
   % extra=halfwidth ("hw" calculated in line 15).
   [FitResults,GOF]=peakfit([x;y],0,0,2,30,hw,10);
   disp('        Peak #      Position     Height     GauWidth       Area       LorWidth')
   disp(FitResults)
   MeasuredAlphas=FitResults(:,4)./FitResults(:,6)