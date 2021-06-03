% Demonstration of a real-time Fourier bandpass filter. It pre-computes
% a simulated signal on line 38. A critical variable in this case is “WindowWidth”
% (line 55), the number of data points taken to compute each filtered
% spectrum. The larger this number, the fewer the number of segments that
% will be generated, but the higher will be the frequency resolution. On a
% standard desktop PC (Intel Core i5 3 Ghz running Windows 10 home), this
% script generates about 35 segments per second with an average data rate
% (points per seconds) of about 34,000 Hz. Smaller segments (i.e. lower
% values of WindowWidth) generate proportionally lower average data rates
% (because the signal stream is interrupted more often to calculate and
% graph the filtered spectrum).
% 
% In this demonstration, the Fourier bandpass filter is used to detect a
% 500 Hz ('f' in line 28) sine wave that occurs in the middle third of a
% very noisy signal (line 32), from about 0.7 sec or 1.3 sec. The filter
% center frequency and width are set in lines 46 and 47.
% Tom O'Haver (toh@umd.edu) 2018
format compact
format short g
figure(1)
clear
clf
% Create simulated signal
SamplingTime=2; % Sampling time in seconds
increment=.00001;
t=0:increment:SamplingTime; % time vector in seconds
SamplingFrequency=length(t)/SamplingTime;
f=500; % signal frequency in cycles per second (Hz)
y=sin(2*pi.*f.*t); % Create sine wave of frequency f
y(1:round(length(t)/3))=0; % Apply sine to middle third of signal
y(round(length(t)/(3/2)):length(t))=0;
noise=2*randn(size(t));
y=y+noise; % add noise
maxn=length(y);
maxy=max(y);
miny=min(y);
plot(t,y)
title('Raw signal')
xlabel('Time, sec.');
ylabel('Amplitude')
pause(2)

 % Set filter variables
CenterFrequency=500; % Center Frequency of bandpass filter, in Hz
FilterWidth=1; % width of bandpass filter, in Hz
WindowWidth=1000; % The number of data points over which each spectrum is performed
PlaySound=1; % Set PlaySound to 1 to play the data as a sound
figure(1)
sn=1;
segment=zeros(1,WindowWidth);
FilteredSignal=[];

% Process the data point-by-point
tic
for n=1:maxn-1 % Display spectra without sound
    ysegment(sn)=y(n);
    tsegment(sn)=t(n);
    sn=sn+1;
    if sn==WindowWidth
        wsegment=ysegment;
        FilteredSegment=FouFilter(wsegment,WindowWidth*increment,CenterFrequency,FilterWidth,5,0);
        FilteredSignal=[FilteredSignal FilteredSegment];
        subplot(2,1,1)
        plot(ysegment)
        subplot(2,1,2)
        plot(FilteredSegment)
        axis([0 WindowWidth miny maxy])
        drawnow
        sn=1;
    end
end
ElapsedTime=toc
SignalDuration=length(y)./SamplingFrequency
TimePerPoint=ElapsedTime/(maxn)
NumFrames=maxn./WindowWidth
FramesPerSecond=NumFrames./ElapsedTime
DataRate=maxn./ElapsedTime
clf
subplot(2,1,1)
plot(t,y)
title('Raw signal Input')
subplot(2,1,2)
plot(t(1:length(FilteredSignal)),FilteredSignal)
xlabel(['Time, sec      Filter frequency = ' num2str(CenterFrequency) '     Filter width = ' num2str(FilterWidth) ] );
ylabel('Amplitude')
title('Real Time Fourier Filter')
grid
