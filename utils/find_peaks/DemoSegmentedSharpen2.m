 % Demonstration of segmented sharpen function.
 % Four identical Gaussian peaks subjected to increasing sharpening
 % using the SegmentedSharpen.m function.
 % Variation using start and end points for sparpening vectors
 % Must have SegmentedSharpen.m and autofindpeaks.m in path.
 format compact
 format short g
 figure(1)
 clf
  x=0:.01:35;
  y=exp(-(x-6).^2)+exp(-(x-14).^2)+exp(-(x-22).^2)+exp(-(x-30).^2+exp(-(x-38).^2));
  NumSegments=5;
  start1=400;
  end1=1400;
  factor1=start1:(end1-start1)/NumSegments:end1; % Increasing values of factor 1
  start2=1.1e6;
  end2=1.4e6;
  factor2=start2:(end2-start2)/NumSegments:end2; % Increasing values of factor 2
  smoothwidth=[3 3 3 3 3]; % Same smooth width for all peaks
  Enhancedsignal=SegmentedSharpen(y,factor1,factor2,smoothwidth);plot(x,y,x,Enhancedsignal,'r')
  plot(x,y,x,Enhancedsignal,'r')
  title('Demonstration of segmented sharpen function.')
  ylabel('y')
  xlabel('Four identical Gaussian peaks (blue) subjected to increasing sharpening (red)')
  OriginalPeaks=autofindpeaks(x,y,30);
  EnhancedPeaks=autofindpeaks(x,Enhancedsignal,30);
  disp('Original Peaks:')
disp('           Peak     Position     Height       Width       Area')
  disp(OriginalPeaks);
  disp(' ')
  disp('Sharpened Peaks:')
disp('           Peak     Position     Height       Width       Area')
  disp(EnhancedPeaks);
  disp(' ')
  disp('% difference between original and sharpened peaks')
disp('           Peak     Position     Height       Width       Area')
  difference=100.*(EnhancedPeaks-OriginalPeaks)./OriginalPeaks;
  disp(difference)

  