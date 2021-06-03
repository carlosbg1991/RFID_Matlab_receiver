 x=1:100;
for w=49:-1:0,
    xrange=(50-w):(50+w);
    y=zeros(size(x));y(xrange)=cos(xrange);
    clf
    subplot(2,1,1)
    plot(x,y)
    xlabel('time')
    ylabel('signal')
    title('The duration of a sine wave pulse is inversly proportional to the width of its frequency spectrum')
    subplot(2,1,2)
    PlotFrequencySpectrum(x',y',1,0,0);
    drawnow
end
pause(.5)
for w=0:1:49,
    xrange=(50-w):(50+w);
    y=zeros(size(x));y(xrange)=cos(xrange);
    subplot(2,1,1)
    plot(x,y)
    xlabel('time')
    ylabel('signal')
    title('The duration of a sine wave pulse is inversly proportional to the width of its frequency spectrum')
    subplot(2,1,2)
    PlotFrequencySpectrum(x',y',1,0,0);
    drawnow
end
pause(.5)
for w=49:-1:0,
    xrange=(50-w):(50+w);
    y=zeros(size(x));y(xrange)=cos(xrange);
    subplot(2,1,1)
    plot(x,y)
    xlabel('time')
    ylabel('signal')
    title('The duration of a sine wave pulse is inversly proportional to the width of its frequency spectrum')
    subplot(2,1,2)
    PlotFrequencySpectrum(x',y',1,0,0);
    drawnow
end
pause(.5)
for w=0:1:49,
    xrange=(50-w):(50+w);
    y=zeros(size(x));y(xrange)=cos(xrange);
    subplot(2,1,1)
    plot(x,y)
    xlabel('time')
    ylabel('signal')
    title('The duration of a sine wave pulse is inversly proportional to the width of its frequency spectrum')
    subplot(2,1,2)
    PlotFrequencySpectrum(x',y',1,0,0);
    drawnow
end