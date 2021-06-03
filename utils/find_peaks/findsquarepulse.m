function S=findsquarepulse(t,y,Threshold)
% function S=findsquarepulse(t,y,Threshold)
% Locates the rectangular pulses in the signal t,y that exceed a y-value of
% "threshold" and determines their start time, average height (relative to
% the adjacent baseline) and width. If the signal is very noisy, some
% preliminary rectangular smoothing (e.g. using fastsmooth.m) before
% calling findsquarepulse.m may be helpful to eliminate false peaks.
% DemoFindsquare.m creates a test signal and calls findsquarepulse.m to
% demonstrate.
% Example: 
% t=1:100;y=sign(sin(t));plot(t,y);S=findsquarepulse(t,y,0)
Pulse=0;
Top=0;
pulsepoints=0;
basepoints=0;
baseline=0;
PulseStart=0;
S=[0 0 0 0];
for k=1:length(t)-1,
    % Points above threshold count as top points
    if y(k) > Threshold,
        Top=Top+y(k);
        pulsepoints=pulsepoints+1;
    end
    % Points below threshold count as base points
    if y(k) < Threshold,
        baseline=baseline+y(k);
        basepoints=basepoints+1;
    end
    % Look for an up transition to signify start of pulse
    if y(k) > Threshold && y(k+1) < Threshold,
        Width=t(k)-PulseStart;
        Height=(Top./pulsepoints)-(baseline./basepoints);
        if Pulse>0,
           S(Pulse,:)=[round(Pulse) PulseStart Height Width];
        end
        baseline=0;
        basepoints=0;
    end
    % Look for a down transition to signify end of pulse
    if y(k) < Threshold && y(k+1) > Threshold,
        Pulse=Pulse+1;
        PulseStart=t(k);
        Top=y(k+1);
        pulsepoints=0;
    end
end