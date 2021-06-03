% Demonstration of the effect of quantization on recorded digitized speech
% and the addition of random white noise before quantization.
warning off all
v=wavread('TestingOneTwoThree.wav');
time=0:1/44001:1.5825;
waveform=v(:,2);

% 256 steps (8 bits), no added noise
NumSteps=256;   
noise=0;    
%
f=NumSteps-1;
s256=(1./f).*round(f.*waveform+noise.*randn(size(waveform)));
figure(1)
plot(time,f*s256)
xlabel('Time, seconds')
ylabel('Digitized waveform')
title('256 steps (8 bits), no added noise: Fair quality')
pause(1)
sound(s256,44000);
wavwrite(s256,44000,'s256');

pause(3)

% 8 steps (3 bits), no added noise
NumSteps=8;  
noise=0;   
%
f=NumSteps-1;
s8=(1./f).*round(f.*waveform+noise.*randn(size(waveform)));
figure(2)
plot(time,f*s8)
xlabel('Time, seconds')
ylabel('Digitized waveform')
title('8 steps (3 bits), no added noise: Poor quality: quantization noise noticable.')
pause(1)
sound(s8,44000);
wavwrite(s8,44000,'s8');

pause(3)


% 2 steps, no added noise
NumSteps=2;  
noise=0;    
%
f=NumSteps-1;
s2=(1./f).*round(f.*waveform+noise.*randn(size(waveform)));
figure(3)
plot(time,f*s2)
xlabel('Time, seconds')
ylabel('Digitized waveform')
title('2 steps, no added noise:  Not intelligible.')
pause(1)
sound(s2,44000);
wavwrite(s2,44000,'s2');

pause(3)


% 2 steps, 0.1 added noise
NumSteps=3;  
noise=0.1;    
%
f=NumSteps-1;
s2n=(1./f).*round(f.*waveform+noise.*randn(size(waveform)));
figure(4)
plot(time,f*s2n)
xlabel('Time, seconds')
ylabel('Digitized waveform')
title('2 steps, 0.1 added noise: Bad quality; barely intelligible.')
pause(1)
sound(s2n,44000);
wavwrite(s2n,44000,'s2n');