 x=1:100;
for w=99.1:-1:.1,
    y=gaussian(x,50,w).*cos(x);
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
for w=.1:1:99.1,
    y=gaussian(x,50,w).*cos(x);
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
for w=99.1:-1:1,
    y=gaussian(x,50,w).*cos(x);
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
for w=1:1:99.1,
    y=gaussian(x,50,w).*cos(x);
    subplot(2,1,1)
    plot(x,y)
    xlabel('time')
    ylabel('signal')
    title('The duration of a sine wave pulse is inversly proportional to the width of its frequency spectrum')
    subplot(2,1,2)
    PlotFrequencySpectrum(x',y',1,0,0);
    drawnow
end