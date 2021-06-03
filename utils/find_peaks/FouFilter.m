function [ry,fy,ffilter,ffy] = FouFilter(y,samplingtime,centerFrequency,filterWidth,shape,mode )
% [ry,fy,ffilter,ffy] = FouFilter( y,samplingtime,centerFrequency,filterWidth,shape,mode )
% Fourier filter function for time-series signal vector y; 'samplingtime'
% is the total duration of sampled signal in sec, millisec, or microsec;
% 'centerfrequency' and 'filterWidth' are the center frequency and width 
% of the filter in Hz, KHz, or MHz, respectively; 'Shape' determines the 
% sharpness of the cut-off. If shape = 1, the filter is Gaussian; as  
% shape increases the filter shape becomes more and more rectangular. 
% Set mode = 0 for band-pass filter, mode = 1 for band-reject (notch) filter.  
% FouFilter returns the filtered signal.
%
% Example: Sine wave in noisy background.
% First half is just noise; sine wave starts halfway through.
% xx=[0:.001:2*pi]';
% signal=sin(20*xx);
% noise=randn(size(xx));
% x=1:2*length(xx)';
% y=[noise;signal+noise]; % sine wave is added halfway through.
% SignalToNoiseRatio=std(signal)/std(noise)
% FilteredSignal=FouFilter(y',1,20,100,5,0);
% subplot(2,1,1)
% plot(x,y);
% title('First half is just noise; sine wave added halfway through')
% subplot(2,1,2)
% plot(x,FilteredSignal);
% title('Signal filtered with FouFilter.m')
%
%  T. C. O'Haver (toh@umd.edu),  version 2, March, 2019
% Correction thanks to Dr. Raphael Attié, NASA/Goddard Space Flight Center

center=centerFrequency*samplingtime; %  center harmonic (fourier component)
width=filterWidth*samplingtime; %  width of filter (in harmonics)

fy=fft(y); % Fourier transform of signal

% lft1=1:(length(fy)/2);
% lft2=(length(fy)/2+1):length(fy);
% 
% ffilter1=ngaussian(lft1,center+1,width,shape);
% ffilter2=ngaussian(lft2,length(fy)-center+1,width,shape);

lft1=0:floor(length(fy)/2)-1;
lft2=floor(length(fy)/2):length(fy)-1;

ffilter1=ngaussian(lft1,center,width,shape);
ffilter2=ngaussian(lft2,length(fy)-center,width,shape);

ffilter=[ffilter1,ffilter2];
if mode==1
    ffilter=1-ffilter; 
end
if length(fy)>length(ffilter), ffilter=[ffilter ffilter(1)];end
ffy=fy.*ffilter;  % Multiply filter by Fourier transform of signal
ry=real(ifft(ffy)); % Recover filter signal from Fourier transform
end

