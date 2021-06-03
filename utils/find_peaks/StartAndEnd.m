% Precision of peak start and stop points by iterative curve fitting.
% Requires peakfit.m, gaussian.m, BuGaussian.m, val2ind.m, and rsd.m in the path.
clear startpoint; clear endpoint;clear centerindex
interval=4; % Change interval to change sampling rate
x=100:interval:450; % define x axis
NumTrials=50; % Number of repeat trials to measure standard deviation
noise=0.010; % random white noise as a fraction of the max peak height
CutOff=0.01; % Fraction of peak height taken as cut-off points
NumPeaks=2; % Number of Gaussian peaks in Model 2.
for n=1:NumTrials
    % Simulate asymmetrical peak with random white noise
    y=BiGaussian(x,240,75,2.5)+noise.*randn(size(x));
    %
    % Select a model to fit the data by un-commenting line 13 or line 15
    % Model 1: Single asymmetrical peak, when shape is known (faster)
    % [FitResults,GOF,baseline,coeff,residuals,xi,yi]=peakfit([x;y],259,500,1,14,2.5,0,0,0,0,0);
    % Model 2: Two overlapping Gaussians (slower but fits wider range of peak shapes)
    [FitResults,GOF,baseline,coeff,residuals,xi,yi]=peakfit([x;y],259,500,NumPeaks,1,0,3.*NumPeaks,0,0,0,0);
    %
    % Calculate start and end points from model, xi, yi, in output arguments
    sizeresults=size(FitResults);
    if sizeresults(1)==1
        sy=yi; % total model for single-peak models
    else
        sy=sum(yi); % total model for multiple-peak models
    end
    maxsy=max(sy); % peak height maximum
    ymin=CutOff.*maxsy; % cut-off amplitude
    % The next 3 statemens use theval2ind(x,val) function, which returns 
    % the index of the element of vector x that is closest to val.
    centerindex(n)=val2ind(sy,maxsy); % point index number at peak center
    centerpoint(n)=xi(val2ind(sy,sy(centerindex(n))));
    startpoint(n)=xi(val2ind(sy(1:centerindex(n)),ymin));
    endpoint(n)=xi(val2ind(sy(centerindex(n):length(xi)),ymin)+centerindex(n));
    %
    % Plot and label figure
    figure(1)
    clf
    plot(x,y,'.',xi,yi,xi,sy,[startpoint(n) startpoint(n)],[0 1.2],[endpoint(n) endpoint(n)],[0 1.2])
    axis([min(x) max(x) -.2 1.2])
    grid
    title('Measurement of peak start and stop points by curve fitting')
    xlabel(['Random noise: ' num2str(100.*noise) '%      CutOff point:  ' num2str(100*CutOff) '%     Start: ' num2str(startpoint(n)) '      End:  ' num2str(endpoint(n)) ])
end
MeanCenter=mean(centerpoint);
PercentRSDcenter=100*rsd(centerpoint);
Meanstart=mean(startpoint);
PercentRSDstart=100*rsd(startpoint);
Meanend=mean(endpoint);
PercentRSDend=100*rsd(endpoint);
disp(' ')
disp('               Mean     % RSD')
disp(['Peak center   '  num2str(MeanCenter) '  ' num2str(PercentRSDcenter) ])
disp(['Start point   '  num2str(Meanstart) '  ' num2str(PercentRSDstart) ])
disp(['End point     '  num2str(Meanend) '  ' num2str(PercentRSDend) ])