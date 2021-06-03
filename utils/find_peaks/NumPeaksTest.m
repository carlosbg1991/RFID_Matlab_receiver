% NumPeaksTest.m demonstrates an attempt to determine the minimum number of
% model peaks needed to fit a set of data, plotting the fitting error vs
% the number of model peaks, and looking for the point at which the fitting
% error reaches a minimum. This script creates a noisy computer-generated
% signal containing a user-selected 3, 4, 5 or 6 underlying peaks, fits to
% a series of models containing 1 to 10 model peaks. The correct number of
% underlying peaks is either the fit with the lowest fitting error, or, if
% two or more fits have about the same fitting error, the fit with the
% least number of peaks. Adjust x, wid, noise, peakshape, and NumPeaks 
% to look like your type of data. Change the shape function in the switch
% statement to match the shape number in line 17.
clear R
x=1:1:100; % (Adjust to look like your type of data)
wid=20; % Width of peaks (Adjust to look like your type of data)
noise=0.02; % Random noise added (Adjust to look like your type of data)
NumTrials=30; % Should be 20 or 30 (larger the better, but it slows down operation)
peakshape=2; % model peak shape: 1=Gaussian, 2=Lorentzian, type 'help peakfit' for others.
NumPeaks=4; % Number of actual underlying peaks in the data
%
switch NumPeaks
    case 3
        y=lorentzian(x,25,wid)+lorentzian(x,40,wid)+lorentzian(x,60,wid);
    case 4
        y=lorentzian(x,25,wid)+lorentzian(x,40,wid)+lorentzian(x,65,wid)+lorentzian(x,80,wid);
    case 5
        y=lorentzian(x,25,wid)+lorentzian(x,40,wid)+lorentzian(x,60,wid)+lorentzian(x,75,wid)+lorentzian(x,90,wid);
    case 6
        y=lorentzian(x,10,wid)+lorentzian(x,25,wid)+lorentzian(x,40,wid)+lorentzian(x,60,wid)+lorentzian(x,75,wid)+lorentzian(x,90,wid);
end
y=y+noise.*randn(size(y));
for TrialPeaks=1:10;
[Results,FitError]=peakfit([x;y],0,0,TrialPeaks,peakshape,0,NumTrials);
 R(TrialPeaks,:)=([TrialPeaks FitError]);
 drawnow
end
clf
 semilogy(R(:,1),R(:,2),'o-'); % plot fitting error vs number of peaks
 xlabel('Number of model peaks')
 ylabel('Percent fitting error')
disp('           Peaks  FitError      R2')
disp(R)
