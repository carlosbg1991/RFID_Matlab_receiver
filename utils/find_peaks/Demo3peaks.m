% A self-contained interactive demonstration of iPeak applied 
% to a data set containing three simple peaks with increasing peak height 
% and peak width. Use this to understand the difference between the 
% variables SlopeThreshold (SlopeT), which discriminates on the basis 
% of peak width, and AmpThreshold (AmpT), which discriminates on the 
% basis of peak amplitude. Peak number and the estimated position, height, 
% and width of each peak is returned in the matrix P. Requires that the
% ipeak.m function be placed in the Matlab path.
% Tom O'Haver (toh@umd.edu). Version 2 August, 2011
warning off all
format compact

% Simulate data set
increment=1;
x=[1:increment:400];
% For each simulated peak, enter the amplitude, position, and width below
amp=[3 4 5];  % Amplitudes of the peaks
pos=[100 200 300];   % Positions of the peaks
wid=[20 40 60];   % Widths of the peaks
Noise=.05;
% A = matrix containing one of the unit-amplitude peak in each of its rows
A = zeros(length(pos),length(x));
for k=1:length(pos)
  if amp(k)>0, A(k,:)=gaussian(x,pos(k),wid(k)); end; % Or you can use any other peak function
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));
y(1)=6;
% Call the interactive findpeaks script
P=ipeak([x;y],0,0.4207,0.0013286,15,37,200,100)

% Increase the AmpThreshold (AmpT) with the A/Z keys to reject smaller
% peaks.  Increase the SlopeThreshold (SlopeT) with the S/X keys to reject 
% broader peaks.  Increase the smooth width (SmoothW) with the D/C keys to 
% reject narrow peaks. Adjust the fit width (FitW) with the F/V keys to
% cover roughly the top half of the peaks.  Press P to display table
% of measured peaks.  Compare to the true peak parameters in lines 17-19.
