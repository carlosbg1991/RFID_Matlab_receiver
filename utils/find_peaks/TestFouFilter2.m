% Demonstration of a Fourier bandpass filter applied to a noisy 
% 100 Hz sine wave signal with the filter center frequency
% swept from 50 to 150 Hz. Requires the FouFilter.m function
% in the Matlab/Octave path.
% Tom O'Haver 2018
SamplingTime=1; % Sampling time in seconds
f=100; % signal frequency in cycles per second (Hz)
noise=1; % Amplitude of white noise added to signal
FilterWidth=10; % width of bandpass filter, in Hz
%
t=0:.0001:SamplingTime; % time vector in seconds
signal=sin(2*pi.*f*t); % Create sine wave of frequency f
whitenoise=noise.*randn(size(t)); % Create random white noise vector
y=signal+whitenoise; % Added random noise throughtout signal
SignalToNoiseRatio=std(signal)/std(noise);
for CenterFrequency=50:150
    FilteredSignal=FouFilter(y,SamplingTime,CenterFrequency,FilterWidth,5,0);
    subplot(2,1,1)
    plot(t,y);
    title(['Noisy Sine Wave      Signal frequency = ' num2str(f) '     SignalToNoiseRatio = ' num2str(SignalToNoiseRatio)])
    subplot(2,1,2)
    plot(t,FilteredSignal);
    axis([0 SamplingTime -1.1 1.1])
    title('Signal filtered with FouFilter.m, swept filter frequency')
    xlabel(['Time, sec      Filter frequency = ' num2str(CenterFrequency) '     Filter width = ' num2str(FilterWidth) ] );
    if CenterFrequency==f,pause(1);end
    drawnow
end
