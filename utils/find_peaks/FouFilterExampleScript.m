% Example: Sine wave in noisy background.
%  First half is just noise; sine wave starts halfway through.
  xx=[0:1:2000*pi]';
  f=0.005; % signal frequency
  signal=sin(f*2*pi*xx);
  noise=randn(size(xx));
  x=1:2*length(xx)';
  y=[noise;signal+noise]; % sine wave is added halfway through.
  samplingtime=1;
  centerfrequency=1/(pi*f);
  frequencywidth=100;
  shape=5;
  mode=0; % for band-pass filter; 1 for band-reject (notch) filter
  FilteredSignal=FouFilter(y',samplingtime,centerfrequency,frequencywidth,shape,mode);
  subplot(2,1,1)
  plot(x,y);
  title('First half is just noise; sine wave starts halfway through')
  subplot(2,1,2)
  plot(x,FilteredSignal);
  title('Signal filtered with FouFilter.m')
  OriginalSignalToNoiseRatio=std(signal)/std(noise);
  disp('------------------------------------------------')
  disp(['Original Signal-To-Noise Ratio = ' num2str(OriginalSignalToNoiseRatio)])
  disp('Signal-to-Noise Ratio of filtered signal:')
  Ly=length(y);
  MeasuredNoise=std(FilteredSignal(1:Ly/2));
  disp(['  Measureed Noise = ' num2str(MeasuredNoise)])
  MeasuredSig=std(FilteredSignal(Ly/2:Ly));
  disp(['  Measureed Signal = ' num2str(MeasuredSig)])
  SNR=MeasuredSig/MeasuredNoise;
  disp(['  Signal-to-noise ratio = ' num2str(SNR)])
   disp(['  Improvement factor = ' num2str(SNR./OriginalSignalToNoiseRatio)])