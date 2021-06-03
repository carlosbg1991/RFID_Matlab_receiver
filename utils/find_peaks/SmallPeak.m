% Measuring a weak peak next to a much larger peak of same basic shape.
% Demonstration of increasing measurement accuracy by position and width
% constraints. Hint: Spread out the figure windows so they can all be seen
% separately and don't overlap. You can change the constants in lines 11 -
% 20 or the peak shape in line 26. Try SmallHeight = .1 and .05; Noise = 0
% and .01; xshift = 0 and 1. Requires gaussian.m, peakfit.m and findpeaks.m
% in the path.
format compact
format short g
clear
increment=.05; % x-axis increment (smaller => more data points)
P1=4; % Peak position of first (larger) interfering peak
P2=5; % Peak position of second (smaller) measured peak
LargeHeight=1; % Height of peak 1 (the uninteresting or interfering peak)
SmallHeight=.05; % Height of peak 2 (the measured peak)
width=1.66; % Full width at half maximum of both peaks
Noise=.01; % Fractional random white noise
xshift=0; % Average unexpected shift in x-axis, causing both peaks to shift
NumTrials=5; % Number of trial fits per signal from different starting points
NumSignals=20; % Number of repeated measurements with independent noise
SignalToNoiseRatio=SmallHeight/Noise % SNR of smaller peak
x=0:increment:10;n=length(x);
for trial=1:NumSignals,
    % Compute signal vector y for this measurement (change shape in line 26)
    shift=xshift.*randn(); % Shift in position of peaks for this measurement
    y=SmallHeight.*gaussian(x,P2+shift,width)+LargeHeight.*gaussian(x,P1+shift,width)+Noise.*randn(size(x));
    
    figure(1)
    F0=peakfit([x;y],0,0,2,1,0,NumTrials,0,0);
    F0=sortrows(F0,3);
    subplot(2,1,1)
    title('Unconstrained iterative fit (variable positions and widths), using peakfit.m: ')
    
    figure(2)
    F1=peakfit([x;y],0,0,2,6,0,NumTrials);
    F1=sortrows(F1,3);
    subplot(2,1,1)
    title('Equal-width iterative fit (variable positions and equal widths), using peakfit.m')
    
    figure(3)
    F2=peakfit([x;y],0,0,2,16,0,NumTrials,0,0,[P1 P2]);
    F2=sortrows(F2,3);
    subplot(2,1,1)
    title('Fixed positions iterative fit (variable widths), using peakfit.m')
 
    figure(4)
    CLSheights=cls(x,y,2,1,[4 5],[width width],0);
    clsp1=LargeHeight.*exp(-(x-P1).^2);
    clsp2=SmallHeight.*exp(-(x-P2).^2);
    clsM=clsp1+clsp2;
    subplot(2,1,1)
    plot(x,y,'b.',x,clsp1,'g',x,clsp2,'g',x,clsM,'r')
    axis([0 10 0 1.2])
    FittingError=100*norm(y-clsM)./(sqrt(n)*max(y));
    xlabel(['Peakheights = '  num2str(CLSheights) '        Fitting error = ' num2str(FittingError) '%' ])
    title('Fixed positions and widths, using CLS.m')
    subplot(2,1,2)
    plot(x,y-clsM,'m.')
    xlabel('Residual plot')
    
    figure(5)
    % Find peak position of this measurement (P(2))
    SlopeThreshold=.003;
    AmpThreshold=.1;
    SmoothWidth=3;
    FitWidth=11;
    Peaktable=findpeaksG(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3); 
    CLSheights2=cls(x,y,2,1,[Peaktable(2) Peaktable(2)+1],[width width],0);
    % Assume peak 1 has this position
    clsp1=LargeHeight.*exp(-(x-Peaktable(2)).^2);
    % Assume peak 2 is also shifted the same amount in this measurement
    clsp2=CLSheights2(2).*exp(-(x-Peaktable(2)-(P2-P1)).^2); 
    clsM=clsp1+clsp2;
    subplot(2,1,1)
    plot(x,y,'b.',x,clsp1,'g',x,clsp2,'g',x,clsM,'r')
    axis([0 10 0 1.2])
    FittingError=100*norm(y-clsM)./(sqrt(n)*max(y));
    xlabel(['Peakheights = '  num2str(CLSheights2) '        Fitting error = ' num2str(FittingError) '%' ])
    title('Position 1 calculated from findpeaks, fixed widths and distance between peaks, using CLS.m')
    subplot(2,1,2)
    plot(x,y-clsM,'m.')
    xlabel('Residual plot')
    
     figure(6)
     start=[Peaktable(2) width Peaktable(2)+(P2-P1) width];
     F5=peakfit([x;y],0,0,2,6,0,NumTrials,start);
     F5=sortrows(F5,3);
     subplot(2,1,1)
     title('Equal-width iterative fit (first guess derived from findpeaks), using peakfit.m')
    
    e0(trial)=100.*(SmallHeight-F0(1,3))./SmallHeight;
    e1(trial)=100.*(SmallHeight-F1(1,3))./SmallHeight;
    e2(trial)=100.*(SmallHeight-F2(1,3))./SmallHeight;
    e3(trial)=100.*(SmallHeight-CLSheights(2))./SmallHeight;
    e4(trial)=100.*(SmallHeight-CLSheights2(2))./SmallHeight;
    e5(trial)=100.*(SmallHeight-F5(1,3))./SmallHeight;
end
Unconstrained=sqrt(norm(e0));
EqualWidths=sqrt(norm(e1));
FixedPositions=sqrt(norm(e2));
FixedPositionsAndWidths=sqrt(norm(e3));
findpeaksP=sqrt(norm(e4));
findpeaksP2=sqrt(norm(e5));
disp('      Unconstr.     EqualW       FixedP      FixedP&W    findpeaksP   findpeaksP2')
disp([Unconstrained EqualWidths FixedPositions FixedPositionsAndWidths findpeaksP findpeaksP2])