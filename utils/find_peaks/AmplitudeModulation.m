% AmplitudeModulation.m is a Matlab/Octave script simulation of modulation
% and synchronous detection, in which the signal created when a light
% beam scans the test sample is modeled as a Gaussian band ('y'), whose
% parameters are defined in the first few lines. As the spectrum of the
% sample is scanned, the light beam is amplitude modulated by the chopper,
% represented as a square wave defined by the bipolar vector 'reference',
% which switches between +1 and -1, shown in the top panel of Figure 1. The
% modulation frequency is many times faster than the rate at which the
% sample is scanned. The light emerging from the sample therefore shows a
% finely chopped Gaussian ('my'), shown in the second panel of Figure 1.
% But the total signal seen by the detector also includes an unstable
% background introduced after the modulation ('omy'), such as lighted
% emitted by the sample itself or detector background, which in this
% simulation this is modeled as a random walk (Appendix O), which seriously
% distorts the signal, shown in the 3rd panel of Figure 1. The detector
% signal is then sent to a lock-in amplifier that is synchronized to the
% reference waveform; the action of the lock-in is to multiply signal by
% the bipolar reference waveform, inverting the signal when the light of
% off and passing it unchanged when the light is on. The result ('dy') is
% shown in the 4th panel of Figure 1. Now this is low-passed filtered to
% remove the modulation frequency, resulting in the recovered signal peak
% 'sdy' shown in the bottom panel of Figure 1.
%   These various signal are compared in figure 2. The underling signal
% peak is shown as the blue line, and the contaminating background is shown
% in black, in this case modeled as a random walk. The total signal that
% would have been seen by the detector before modulation was used is shown
% in green; the signal distortion is evident, and any attempt to measure
% the signal peak in that signal would be greatly in error. The signal
% recovered by the modulation and lock-in system is shown in red and
% overlayed with the original signal peak in blue for comparison. The
% script also uses the peakfit.m function to measure the peak parameters in
% the original unmodulated total signal (green line) and in the modulated
% recovered signal and to compute the relative percent errors by both
% methods and print them out in a table, Requires fastsmooth.m and
% peakfit.m in the path.(http://tinyurl.com/cey8rwh). Version 2, June 2018,
% modified for Matlab version R2017. Tom O'haver, 2018

% Changable parameters:
PeakHeight=1; % Signal peak height
PeakPosition=100; % Signal peak position
PeakWidth=40; % Signal peak width
SmoothWidth=60; % post-detection smooti,g filter
NoiseLevel=0.01; % random white noise in signal
Frequency=6; % Modulation frequency
phaseshift=0; % Phase shift between modulaiton and reference waveform (ideally zero)

% Create signal and noisy baseline
x=1:.1:200;
noise=NoiseLevel*randn(size(x));
baseline=10.*noise+cumsum(noise); % random walk = cumulative sum of white noise
y=PeakHeight.*gaussian(x,PeakPosition,PeakWidth);
SignalToNoiseRatio=PeakHeight/std(baseline)
oy=y+baseline; % Observed signal on baseline

% Single Gaussian peak fit, baseline mode 1
[FitResults1,GOF1]=peakfit([x;oy],0,0,1,1,0,10,0,1,0,0,1);

% Modulation 
reference=sign(sin(Frequency.*x)); % 50 percent duty-cycle square wave modulation
my=y.*(1+sign(sin(Frequency.*x+phaseshift))); % Modulated signal
omy=my+baseline;     % Modulated signal on (unmodulated) baseline

% Synchronous Detection (e.g. Lock-in amplifier)
dy=omy.*reference;
sdy=fastsmooth(dy,SmoothWidth,2,1); % Lowpass filter ro remove modulation frequency

% Measure peak parameters: single Gaussian peak fit, baseline mode 1
[FitResults2,GOF2]=peakfit([x;sdy],0,0,1,1,0,10,0,1,0,0,1);

figure(1)
clf
subplot(5,1,1);plot(x,reference);
title('Reference waveform')
subplot(5,1,2);plot(x,my);
title('Modulated peak')
subplot(5,1,3);plot(x,omy);
title('Observed modulated signal, with baseline')
subplot(5,1,4);plot(x,dy)
title('Output of synchronous detector (lock-in amplifier)')
subplot(5,1,5);plot(x,sdy)
title('Smoothed synchronous detector output, to remove modulation frequency')

figure(2)
clf
plot(x,y,'bl',x,baseline,'k',x,y+baseline,'g',x,sdy,'r')
title('Blue = underlying signal     Black = baseline     Green = signal+baseline    Red = recovered signal')

% Tabulate relative percent errors of peak position, height, and width
PositionError1=abs(100.*(PeakPosition-FitResults1(1,2))./PeakPosition);
HeightError1=abs(100.*(PeakHeight-FitResults1(1,3))./PeakHeight);
WidthError1=abs(100.*(PeakWidth-FitResults1(1,4))./PeakWidth);
PositionError2=abs(100.*(PeakPosition-FitResults2(1,2))./PeakPosition);
HeightError2=abs(100.*(PeakHeight-FitResults2(1,3))./PeakHeight);
WidthError2=abs(100.*(PeakWidth-FitResults2(1,4))./PeakWidth);
disp('       Position Error  Height Error  Width Error')
disp(['Original:  ' num2str([PositionError1 HeightError1 WidthError1] )] )
disp(['Modulated: ' num2str([PositionError2 HeightError2 WidthError2] )] )