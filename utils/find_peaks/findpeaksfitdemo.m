% Demonstration of findpeaksfit for finding and fitting the peaks in 150
% signals, each of which have 1 to 3 noisy Lorentzian peaks in variable
% locations.
pausetime=0;
peakwidth=10;
x=1:.2:100;
for k=1:50,
    y=lorentzian(x,10+30.*rand(),peakwidth)+.05.*randn(size(x));
    if rand()>.5,
        [findpeaksr,peakfitr]=findpeaksfit(x,y,8e-005,0.5,17,31,3,2,0,1,0,0,1);
        drawnow
        pause(pausetime)
    end
    y=lorentzian(x,10+30.*rand(),peakwidth)+lorentzian(x,50,peakwidth)+.05.*randn(size(x));
    if rand()>.5,
        [findpeaksr,peakfitr]=findpeaksfit(x,y,8e-005,0.5,17,31,3,2,0,1,0,0,1);
        drawnow
        pause(pausetime)
    end
    y=lorentzian(x,10+30.*rand(),peakwidth)+lorentzian(x,50,peakwidth)+lorentzian(x,60+30.*rand(),peakwidth)+.05.*randn(size(x));
    if rand()>.5,
        [findpeaksr,peakfitr]=findpeaksfit(x,y,8e-005,0.5,17,31,3,2,0,1,0,0,1);
        drawnow
        pause(pausetime)
    end
end