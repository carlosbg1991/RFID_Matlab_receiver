% Demonstration of a Fourier bandpass filter applied to a noisy 100 Hz sine
% wave which appears only in the middle third of the signal record.
% Requires the FouFilter.m function in the Matlab/Octave path.
% Try reducing the filter width to see how the transient response changes.
% Tom O'Haver 2018
SamplingTime=1; % Sampling time in seconds
f=100; % signal frequency in cycles per second (Hz)
noise=2; % Amplitude of white noise added to signal
CenterFrequency=100; % Center Frequency of bandpass filter, in Hz
FilterWidth=20; % width of bandpass filter, in Hz
% Create signal
t=[0:.0001:SamplingTime]; % time vector in seconds
y=sin(2*pi.*f*t); % Create sine wave of frequency f
y(1:round(length(t)/3))=0; % Apply sine to middle third of signal
y(round(length(t)/(3/2)):length(t))=0;
y=y+noise.*randn(size(t)); % Added random noise throughtout signal
% Apply Fourier filter and plot results
FilteredSignal=FouFilter(y,SamplingTime,CenterFrequency,FilterWidth,5,0);
subplot(2,1,1)
plot(t,y);
title('Sine wave appears in middle third of signal.')
subplot(2,1,2)
plot(t,FilteredSignal);
axis([0 SamplingTime -1 1])
title('Signal filtered with FouFilter.m')
xlabel(['Time, sec      Filter frequency = ' num2str(CenterFrequency) '     Filter width = ' num2str(FilterWidth) ] );

