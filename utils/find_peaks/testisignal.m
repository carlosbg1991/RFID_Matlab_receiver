x=[0:.005:2];y=humps(x);Data=[x;y];
isignal(Data);
subplot(2,1,1);title('Example of Humps function.') 

pause(1)
isignal(Data,0.9,0.5);
subplot(2,1,1);title('Zoom to the peak at x=0.9') 

pause(1)
isignal(Data,0.9,0.5,1,3,1,4);
subplot(2,1,1);title(' 4th derivative of the peak at x=0.9') 

pause(1)
isignal(Data,0.3,0.5,1,3,1,0,1,220,5400);
subplot(2,1,1);title('Peak sharpening applied to the peak at x=0.3.') 

pause(1)
x=[0:.01:20];
y=exp(-(x-4).^2)+exp(-(x-9).^2)+exp(-(x-13).^2)+exp(-(x-15).^2);
isignal([x;y]);
subplot(2,1,1);title('Peak area measured by perpendiclar drop method ') 
pause(1)
isignal([x;y],4,5.91);
subplot(2,1,1);title('First peak. ') 
pause(1)
isignal([x;y],8.7,4.8);
subplot(2,1,1);title('Second peak.') 
pause(1)
isignal([x;y],12.53,2.9);
subplot(2,1,1);title('Third peak.') 
pause(1)
isignal([x;y],16.42,4.94);
subplot(2,1,1);title('Fourth peak.') 
pause(1)

x=-5:.01:5;
y=exp(-(x).^2);for n=1:1000,if randn()>2,y(n)=rand()+y(n);end,end;
isignal([x;y],0,10);
subplot(2,1,1);title('Single peak with random spikes.') 
pause(1)
isignal([x;y],0,10,0,0,0,0,0,0,0,0,3,0);
subplot(2,1,1);title('Median filter (M key) to remove spikes.') 

pause(1)
x=0:.1:60; y=sin(x)+sin(10.*x);
[pY,SpectrumOut]=isignal([x;y],30,30,4,3,1,0,0,1,0,0,0,1);
clf;
subplot(2,1,1);
plot(x,y);
subplot(2,1,2);
plot(SpectrumOut) 
xlabel('Frequency Spectrum')
subplot(2,1,1);
title('Direct entry into frequency spectrum mode.');


disp('End of test')