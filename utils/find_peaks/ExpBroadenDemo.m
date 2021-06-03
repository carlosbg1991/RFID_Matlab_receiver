   clf
x=[0:.01:18];
   y=exp(-(x-4).^2)+exp(-(x-9).^2)+exp(-(x-12).^2)+exp(-(x-13.7).^2);
   t=50; % time constant (in data points)
   BroadenedSignal=ExpBroaden(y',-t);
   plot(x,y,x,BroadenedSignal,'r')
   xlabel('Blue: original signal.    Red: broadened signal.')