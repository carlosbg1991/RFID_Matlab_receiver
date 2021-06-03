% The script RealTimeFrequencySpectrumWindow.m computes and plots the
% Fourier frequency spectrum of a signal. Like the scripts above, it loads
% the simulated real-time data from a “.mat file” (in line 31) and then
% accesses that data point-by-point in the processing loop (lines 48-74). A
% critical variable in this case is “WindowWidth” (line 37), the number of
% data points taken to compute each frequency spectrum. The larger this
% number, the fewer the number of spectra that will be generated, but the
% higher will be the frequency resolution. On a standard desktop PC (Intel
% Core i5 3 Ghz running Windows 10 home), this script generates about 50
% spectra per second with an average data rate (points per seconds) of
% about 50,000 Hz. Smaller spectra (i.e. lower values of WindowWidth)
% generate proportionally lower average data rates (because the signal
% stream is interrupted more often to calculate and graph a spectrum).
% 
% If the data stream is an audio signal, it's also possible to play the
% sound through the computer's sound system synchronized with the display
% of the frequency spectra; to do this, set PlaySound=1 in line 38. Each
% segment of the signal is played, in line 33, just before the spectrum of
% that segment is displayed. The sound reproduction will not be not
% perfect, because of the very slight delay while the computer computes and
% displays the spectrum before going on to the next segment. You can adjust
% the pitch of the sound line 40 to make the sound the correct pitch. In
% this script, the data file is in fact an audio recording of a short
% segment of the 'Hallelujah Chorus' from Handel's Messiah, which is the
% source of the figure above (with a WindowWidth of 1024).
% Tom O'Haver (toh@umd.edu) 2018
format compact
format short g
figure(1)
clf
load handel.mat % Contains signal vector (y) and sampling frequency
SamplingFrequency=Fs; % Enter sampling frequency here.
maxf=SamplingFrequency./2; % Maximum (Nyquist) frequency
maxn=73113; %  maximum number of data points acquired before stopping.
maxy=10000; % Maximum amplitude (y axis) on graph
startn=1; % Number of data points skipped in the begining
WindowWidth=1024; % Length of segments for which spectra are calculated
PlaySound=1; % Set PlaySound to 1 to play the data as a sound
time=(1:WindowWidth)./SamplingFrequency;
figure(1)
axis([0 maxf 0 maxy]);
xlabel('Frequency, Hz')
ylabel('Amplitude')
title('Real Time Frequency Spectrum')
sn=1;
segment=zeros(1,WindowWidth);
tic
if PlaySound
    for n=startn:maxn-1 % Play sound with spectra
        segment(sn)=y(n);
        sn=sn+1;
        if sn==WindowWidth
            sound(segment,SamplingFrequency)
            SpectrumMatrix=FrequencySpectrum(time,segment);
            plot(SpectrumMatrix(1,:),SpectrumMatrix(2,:))
            axis([0 maxf 0 maxy]);
            drawnow
            sn=1;
        end
        pause(.093./WindowWidth) % Adjust number to control pitch of sound.
    end
else
    for n=startn:maxn-1 % Display spectra without sound
        segment(sn)=y(n);
        sn=sn+1;
        if sn==WindowWidth
            SpectrumMatrix=FrequencySpectrum(time,segment);
            plot(SpectrumMatrix(1,:),SpectrumMatrix(2,:))
            axis([0 maxf 0 maxy]);
            drawnow
            sn=1;
        end
    end
end
ElapsedTime=toc
SignalDuration=length(y)./SamplingFrequency
TimePerPoint=ElapsedTime/(maxn-startn)
NumSpectra=maxn./WindowWidth
SpectraPerSecond=NumSpectra./ElapsedTime
DataRate=maxn./ElapsedTime
xlabel('Frequency, Hz')
ylabel('Amplitude')
title('Real Time Frequency Spectrum')
grid
