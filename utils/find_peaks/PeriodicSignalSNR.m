% Matlab/Octave script demonstrating the estimation of the peak-to-peak and
% root-mean-square signal amplitude and the signal-to-noise ratio of a
% periodic waveform, estimating the noise by looking at the time periods
% where its envelope drops below a threshold. Needs fastsmooth.m and 
% PlotFrequencySpectrum in the path.
format compact
% Load the wav file into a Matlab variable
v=wavread('TestingOneTwoThree.wav');
time=0:1/44001:1.5825; % 44001 is the sampling rate; 1.5825 is the duration of the signal
waveform=v(:,1)'+v(:,2)';
% Set the threshold below which is mostly noise.
threshold=.03

% load george-dataShort.csv
% waveform=george_dataShort;
% slength=length(waveform);
% time=(1:slength)./8000;

disp(' ')
% Compute the envelope by smoothing the absolute value of the waveform
SmoothWidth=500;
SmoothType=2;
envelope=fastsmooth(abs(waveform),SmoothWidth,SmoothType); 
bz=zeros(size(time));
j=1;
for i=1:length(time),
    if envelope(i)<threshold,
        bz(i)=waveform(i); % bz is signal with segments below threshold zeroed out
        cbs(j)=bz(i); % cbs is concatenated background segments
        j=j+1;
    end;
end;
figure(1)
clf
plot(time,waveform,'c',time,envelope,'r')
xlabel('time, seconds')
ylabel('Amplitude')
title('Blue: Signal waveform      Red: Envelope')
figure(2)
plot(time,waveform,'c',time,bz,'r')
xlabel('time, seconds')
ylabel('Signal waveform')
title(['Blue: Signal waveform    Red: portions falling below the threshold of ' num2str(threshold) ' units' ])
PeakToPeak=max(waveform)-min(waveform); % Peak to peak signal amplitude
rms=sqrt(mean(waveform.^2)); % Root mean square (RMS) signal amplitude
noise=std(cbs); % Standard deviation of extracted noise segments
PeakToPeak_SignalToNoiseRatio=PeakToPeak./noise
RMS_SignalToNoiseRatio=rms./noise
% Inspect power spectrum of concatenated noise segments
figure(3)
clf
% xn = time axis of concatenated background segments
xn=[0:time(2)-time(1):(length(cbs)-1)*(time(2)-time(1))]; 
% cbs is waveform of concatenated background segments
plot(xn,cbs)
xlabel('time, seconds')
ylabel('amplitude')
title('Plot of  concatenated background segments containing mostly noise.')
figure(4)
PowerSpectrum=PlotFrequencySpectrum(xn',cbs',0,0,1);
plot(PowerSpectrum(1:25,1),PowerSpectrum(1:25,2))
xlabel('Frequency, Hz')
ylabel('amplitude')
title('Power spectrum of concatenated noise segments')
% ifilter(time,waveform,1266.7064,700,10,2,'Band-pass');
% isignal(time,waveform);
% Peaks=findpeaksx(time,abs(waveform),1e-6,.05,200,15);figure(2);plot(Peaks(:,2),Peaks(:,3))