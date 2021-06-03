function PSM=PlotSegFreqSpect(x,y,NumSegments,MaxHarmonic,logmode)
% function PSM=PlotSegFreqSpect(x,y,NumSegments,MaxHarmonic,logmode)
% Breaks y into 'NumSegments' equal-length segments, computes the power
% spectrum of each segment, and plots the result of the first 'MaxHarmonic'
% Fourier components as a contour plot. If the number of segments and of
% data points is such that the last segment is incomplete, it is discarded.
% Returns the power spectrum matrix (time-frequency-amplitude) as a matrix
% PSM of size NumSegments x MaxHarmonic. If logmode=1, it computes and
% plots the base10 log of the amplitudes. Figure window 2 shows the 3-D
% mesh plot of the power spectrum matrix, whcih can be rotated by dragging
% the pointer.
% Version 1. (c) Tom O'Haver December, 2019.
%
% EXAMPLE 1: Segmented power spectrum of passing automobile horn
% See https://terpconnect.umd.edu/~toh/spectrum/CaseStudies.html#Harmonic
% load horngoby
% PSM=PlotSegFreqSpect(t,doppler,40,350,0); % 40 segments, 240 harmonics
%
% EXAMPLE 2: Segmented power spectrum of spoken phrase "Testing 1, 2, 3"
% load testing123;
% t=0:1/44001:1.5825; 
% y=testing123(:,1); % Select one of the two stereo channels.
% PSM=PlotSegFreqSpect(t,y,40,160,0); % 40 segments, 160 harmonics.
%
% EXAMPLE 3: Two Gaussians with constant sinewave noise, displayed in log
% mode.
% x=1:12000;
% y=90*gaussian(x,3000,500)+90*gaussian(x,7500,500);
% y=y+sin(x/2); % Add sinewave "noise"
% clf
% PSM=PlotSegFreqSpect(x,y,40,30,1);
%
% Related functions: FouFilter.m, PlotSegFrequencySpectrum.m, 
% SegmentedFouFilter.m,  isignal.m.

sizex=size(x);if sizex(1)>sizex(2);x=x';end
sizey=size(y);if sizey(1)>sizey(2);y=y';end
if NumSegments<2;disp('Number of segments must be >2');NumSegments=2;end
if MaxHarmonic<1;disp('Number of harmonics must be >1');end
PSM=[];
ly=length(y);
SegLength=round(ly./NumSegments);
fy=zeros(1,ly);
sy=zeros(1,ly);
for Segment=2:NumSegments
    startindex=(1+(Segment-1)*SegLength);
    endindix=startindex+SegLength-1;
    if endindix>ly,endindix=ly;end
    indexrange=startindex:endindix;
    fy(indexrange)=fft(y(indexrange));
    sy(indexrange)=fy(indexrange) .* conj(fy(indexrange));
    try
      PSM=[PSM;sy(indexrange+1)];
    catch 
      NumSegments=NumSegments-1;    
    end
end
if logmode
   PSM=log10(PSM);
end
figure(1)

xlabel('Time Segments')
ylabel('Frequency (Harmonics)')
if logmode
   title('Segmented Power Spectrum, log z-axis amplitude scale')
else
   title('Segmented Power Spectrum, linear z-axis amplitude scale')
end
grid
figure(2);mesh(PSM(1:NumSegments,1:MaxHarmonic)'); 