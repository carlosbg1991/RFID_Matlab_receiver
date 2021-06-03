% Script to demonstrate the peak symmetrization mode of iSignal (Shift-Y)
format compact
format short g
x=[19:1/80:60]'; % time
a0= 0.6; % peak 1
a1= 25; % retention time
a2=1; % width of the peak
a3=2; % distortion
EMG1=(a0/2*a3).*exp((a2.^2/(2*a3.^2))+(a1-x)/a3).*(erf((x-a1)/(sqrt(2)*a2)-a2/(sqrt(2)*a3))+(a3/abs(a3))); 
a0= 0.3; %  peak 2
a1= 32; % retention time
a2=1.5; % width of the peak
a3=2; % distortion
EMG2=(a0/2*a3).*exp((a2.^2/(2*a3.^2))+(a1-x)/a3).*(erf((x-a1)/(sqrt(2)*a2)-a2/(sqrt(2)*a3))+(a3/abs(a3)));
clf
TrueArea1=trapz(x,EMG1)
TrueArea2=trapz(x,EMG2)
y=EMG1+EMG2;
isignal(x,y);
disp('Press Shift-Y, enter "1", and press Enter. Then press 1 and/or 2 keys')
disp('to adjust the symmetrization factor until the baseline after the peak')
disp('is a low as possible without going negative. Press J to measured peak')
disp('areas, widths, and symmetry before and after. ')
disp(' ')
disp('Peak areas by perpendicular drop using autopeaks.m')
autopeaksresults=autopeaks(x,y,3);
PDArea1=autopeaksresults(1,5)
PDArea2=autopeaksresults(2,5)