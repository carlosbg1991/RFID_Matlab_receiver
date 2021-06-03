% testmeasurepeaks script: Runs all the examples in the measurepeaks help
% file, with a 1-second pause between each one.
format short g
format compact
disp(' ')
disp('Test of measurepeaks.m')
disp('Example 1: sin(x).^2 has theoretical peaks at x=0.5pi, 1.5pi, 2.5pi,')
disp('3.5pi..., with peak heights of 1.0 and peak areas of pi/2 = 1.5708.')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
x=[0:.01:20]';y=sin(x).^2;
disp(measurepeaks(x,y,0,0,5,5,1))
title('Example 1: sin(x).^2 has peak heights of 1.0 and peak areas of pi/2 = 1.5708.')
pause(1)
disp(' ')

disp('Example 2: The built-in "humps" function has two peaks, no noise.')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(measurepeaks(0:.01:2,humps(0:.01:2),0,0,1,1,1))
title('Example 2: The built-in "humps" function has two peaks, no noise.')
pause(1)
disp(' ')

disp('Example 3: Series of peaks that get progressively taller and wider.')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
x=[0:.01:5]';y=x.*sin(x.^2).^2;
disp(measurepeaks(x,y,0,0,5,5,1))
title('Example 3: Series of peaks that get progressively taller and wider.')
pause(1)
disp(' ')

disp('Example 4: Like example 3, with random white noise added.')
x=[0:.01:5]';y=.1.*randn(size(x))+x.*sin(-x.^2).^2;
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(measurepeaks(x,y,.001,.5,15,15,1))
title('Example 4: Like example 3, with random white noise added.')
pause(1)
disp(' ')

disp('Example 5; Like example 3, with added rising baseline.')
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
x=[0:.01:5]';y=x+x.*sin(x.^2).^2;
disp(measurepeaks(x,y,0,0,5,5,1))
title('Example 5; Like example 3, with added rising baseline.')
pause(1)
disp(' ')

disp('Example 6: Gaussian on linear baseline, theoretical area = 1.7725')
x=[0:.1:10];y=2+x/10+exp(-(x-5).^2)+.01.*randn(size(x));
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(measurepeaks(x,y,0.0001,0.3,3,5,1))
title('Example 6: Gaussian on linear baseline, theoretical area = 1.7725')
xlabel('Peak-valley and Tan skim are most accurate height and area measures.')
pause(1)
disp(' ')

disp('Example 7: Two overlapping Gaussians, zero baseline, theoretical area = 1.7725')
x=[0:.05:10];y=exp(-(x-6).^2)+exp(-(x-3.5).^2)+.01.*randn(size(x));
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(measurepeaks(x,y,0.0001,0.5,8,6,1))
title('Example 7: Two overlapping Gaussians, zero baseline, theoretical area = 1.7725')
disp('PeakMax and Perp drop are most accurate height and area measures.')
xlabel('PeakMax and Perp drop are most accurate height and area measures.')
pause(1)
disp(' ')

disp('Example 8: Narrow Gaussian peak on sloping linear baseline, followed by')
disp('much broader peak, theoretical heights 1 and 1, areas 1.59 and 31.9.')
x=1:.2:150;
y=(150-x)./200+gaussian(x,30,1.5)+gaussian(x,90,30);
disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
disp(measurepeaks(x,y,[0.001 .00005],[.6 .6],[5 120],[8 150],1))
title('Example 8: Two Gaussians on sloping baseline, theoretical heights 1 and 1, areas 1.59 and 31.9.')

