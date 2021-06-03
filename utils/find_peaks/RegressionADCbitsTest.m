% Demonstration of the effect of analog-to-digital converter resolution
% (defined in line 9) on Classical Least Squares (multilinear regression)
% of two very closely-spaced "noiseless" overlapping Gaussian peaks.
% Normally, the random noise (line 10) produces a far greater error.
format compact
format short g
clear
clf
ADCbits=8; % <<<<< CHANGE ADC resoltion (number of bits) here.
noise=.001; % Also try different anounts of random noise here
t=[1:10:1000];
amp=[.5 .01];  % Amplitudes of the peaks (two different values to test dynamic range)
pos=[500 550];   % Positions of the peaks (closer = more difficult)
wid=[200 200];   % Widths of the peaks (wider = more difficult)
for trial=1:20, % Multiple trials to get better estimate of average error
    % A = calibration matrix containing one unit-amplidude peak in each of its rows
    A = zeros(length(pos),length(t));
    for k=1:length(pos)
        A(k,:)=round((2.^ADCbits).*gaussian(t,pos(k),wid(k)))./(2.^ADCbits); % You can use any peak function here
    end
    spectrum=amp*A;
    clf
    % Simulation
    ObservedSpectrum = round((2.^ADCbits).*spectrum + noise.*randn(size(t)))./(2.^ADCbits); % adds noise
    hold on
    % Regress ObservedSpectrum onto calibration matrix of reference spectra "A".
    MeasuredAmp = ObservedSpectrum*A'*inv(A*A'); % "A" is the calibration matrix
    figure(1);plot(t,ObservedSpectrum,'m.')
    plot(t,amp(1).*A(1,:),t,amp(2).*A(2,:))
    plot(t,MeasuredAmp(1).*A(1,:),t,MeasuredAmp(2).*A(2,:))
    plot(t,MeasuredAmp(1).*A(1,:)+MeasuredAmp(2).*A(2,:),'k')
    title(['Analysis of mixture of peaks by multiple regression, using ' num2str(ADCbits) ' bit analog-to-digital converter ' ])
    xlabel('Dots: Measured data    Blue = Component 1    Green = Component 2    Black = Mixture of Components 1 and 2')
    hold off
    MeasuredAmp=MeasuredAmp;
    Error = MeasuredAmp-amp;
    RelativePercentError(trial,:) = abs(100.*Error./MeasuredAmp);
end % for trial
    disp(' ')
SeparationToWidthRatio=(pos(2)-pos(1))./mean(wid)
RandomNoise=noise
disp('   Component 1   Component 2')
MeasuredAmp=MeasuredAmp
RelativePercentError=mean(RelativePercentError)