% Script to create and plot a computer-generated signal with multiple peaks
% on a variable baseline plus variable random noise. You may change the
% lines here marked by <<< to modify the character of the signal peaks,
% baseline, and noise. For example, you can change the peak shape in line
% 17, peak positions, heights, and widths in lines 12-14, the sampling
% density in line 10, baseline shape in line 30, baseline amplitude in line
% 18, noise type in line 37, and noise amplitude in line 16.
%  Tom O'Haver, toh@umd.edu    2015
%
increment=.4; %  <<< Change data density here (total number of points)
x=1:increment:100; % Range of x values
Heights=[5 4 3 2 1]; % <<< Height of each peak above the baseline
Positions=[15 30 50 70 85]; %  <<< Peak positions (between 1 and 100)
Widths=[6 6 6 6 6]; % <<< Width (FWHM) of each peak
Extras=[0 0 0 0 0]; %  <<< Needed only for variable shapes (Voigt, etc)
noise=.1;  %  <<< Change noise amplitude here;  
peakshape=1;  %  <<< Change peak shape here: 1=Gaussian, 2=Lorentzaian, etc.
BaselineIntensity=2; %  <<< Change baseline amplitude here (0-10)
extra=0; % Needed only for variable shapes (Voigt, etc)
NumTrials=5; % <<< Normally 1 - 10
width=mean(Widths);
NumPeaks=length(Heights);
% Baseline functions
FlatBaseline=ones(size(x));
LinearBaseline=(1-x/max(x));
ExpBaseline=exp(-x./(10*width));
GaussBaseline=gaussian(x,0,20*width);
% Baseline may be FlatBaseline, LinearBaseline, ExpBaseline, GaussBaseline
Baseline=BaselineIntensity.*ExpBaseline; %  <<< Change baseline type here
y=modelpeaks(x,NumPeaks,peakshape,Heights,Positions,Widths,Extras)+Baseline;
WhiteNoiseArray=randn(size(x));
PinkNoiseArray=pinknoise(length(x));
PropNoiseArray=propnoise(y);
clf
% Add noise to signal. noise may be white, pink, blue, proportional, or
% sqrt.
y=y+noise.*PropNoiseArray; % <<< Change noise type here (WhiteNoiseArray,PinkNoiseArray, or PropNoiseArray
plot(x,y)