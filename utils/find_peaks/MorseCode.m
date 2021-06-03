% Demonstration of Fourier narrow bandpass filter. MorseCode.m uses iFilter
% to demonstrate the abilities and limitations of Fourier filtering. It
% creates a pulsed fixed frequency (0.05) sine wave that spells out “SOS”
% in Morse code (dit-dit-dit/dah-dah-dah/dit-dit-dit), adds random white
% noise so that the SNR is very poor (about 0.1 in this example). The white
% noise has a frequency spectrum that is spread out over the entire range
% of frequencies; the signal itself is concentrated mostly at a fixed
% frequency (0.05) but the presence of the Morse Code pulses spreads out
% its spectrum over a narrow frequency range of about 0.0004. This suggests
% that a Fourier bandpass filter tuned to the signal frequency might be
% able to isolate the signal from the noise. As the bandwidth is reduced,
% the signal-to-noise ratio improves and the signal begins to emerges from
% the noise until it becomes clear, but if the bandwidth is too narrow, the
% step response time is too slow to give distinct “dits” and “dahs”. The
% step response time is inversely proportional to the bandwidth. (Use the ?
% and " keys to adjust the bandwidth. Press 'P' or the Spacebar to hear the
% sound). Watch video on YouTube at https://youtu.be/agjs1-mNkmY
clear
clf
t=1:1:4000;
f=.05; % Signal frequency
noise=4; % RMS random noise
s=sin(2*pi*f*t); 
space=zeros(size(t)); % Small interval of silence
dit=[space s ];
dash=[space s s s s s];
ess=[space dit dit dit space]; % the letter "S" in Morse Code
oh=[space dash dash dash space];  % the letter "O" in Morse Code
sos=[space ess oh ess space];  % "SOS" in Morse Code, surrounded by spaces
signal=std(sos);
SNR=signal/std(noise.*randn(size(sos))) % Signal-To-Noise ratio
nsos=sos+noise.*randn(size(sos));  % Add lots of random white noise 
nsos=nsos./max(sos);
disp(' ')
disp('Look at the title of the figure')
fsos=ifilter(1:length(nsos),nsos,0.05,2,9,1,'Band-pass');
subplot(3,1,1);title('Bandwidth = 2. Bandwidth is too wide.')
sound(fsos,44000);
pause(1)
fsos=ifilter(1:length(nsos),nsos,0.05,.1,9,1,'Band-pass');
subplot(3,1,1);title('Bandwidth = .1. Slightly better.')
sound(fsos,44000);
pause(1)
fsos=ifilter(1:length(nsos),nsos,0.049906,0.005,9,1,'Band-pass');
subplot(3,1,1);title('Bandwidth = 0.005. Even better.')
sound(fsos,44000);
pause(1)
fsos=ifilter(1:length(nsos),nsos,0.05,0.0008,9,1,'Band-pass');
subplot(3,1,1);title('Bandwidth = 0.0008. Much better')
sound(fsos,44000);
pause(1)
fsos=ifilter(1:length(nsos),nsos,0.05,0.00012,9,1,'Band-pass');
subplot(3,1,1);title('Bandwidth = 0.00012. Bandwith is too narrow: response time is too slow.')
sound(fsos,44000);
