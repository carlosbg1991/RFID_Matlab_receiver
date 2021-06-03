x=50:150;
 for Trial=1:10,
     PeakHeight=5.*Trial-1;
     y=PeakHeight.*gaussian(x,100,100)+10.*randn(size(x));
     [H,P,W]=gaussfit(x,y);
     Height(Trial)=H;Position(Trial)=P; Width(Trial)=W; 
 end
 plotit(0:9,Height,1);