function R=testnumpeaks(x,y,peakshape,extra,NumTrials,MaxPeaks)
% function R=testnumpeaks(x,y,peakshape,extra,NumTrials,MaxPeaks)
% Simple test to estimate the number of model peaks required to fit an x,y
% data set. Fits data x,y, with shape "peakshape", with optional extra
% shape factor "extra", with NumTrials repeat per fit, up to a maximum of
% "MaxPeaks" model peaks, displays each fit and graphs fitting error vs
% number of model peaks. If two or more numbers of peaks give about the
% same error, take the smaller number.  toh@umd.edu Feb. 27, 2014, revised
% January, 2108
%
% Example:
% x=0:.01:10;
% y=exp(-(x-4).^2) + exp(-(x-5.3).^2) + exp(-(x-6.5).^2)+.05.*randn(size(x));
% testnumpeaks(x,y,1,0,5,5)
%
for TrialPeaks=1:MaxPeaks
 [~,FitError]=peakfit([x;y],0,0,TrialPeaks,peakshape,extra,NumTrials);
 R(TrialPeaks,:)=([TrialPeaks FitError]);
 drawnow
end
clf
 plot(R(:,1),R(:,2),'o-'); % plot fitting error vs number of peaks
 xlabel('Number of model peaks')
 ylabel('Percent fitting error')
disp('           Peaks  FitError      R2')
disp(R)