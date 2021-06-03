function ipower
% Keyboard-controlled interactive power spectrum demonstrator, useful for 
% teaching and learning about the power spectra of different types of
% signals and the effect of signal duration and sampling rate. Keystrokes
% allow you to select the type of signal, the total duration of the 
% signal, the sampling rate, and the global variables f1 and f2 which are 
% used in different ways in the different signals.  When the Enter key is 
% pressed, the signal (y) is sent to the Windows WAVE audio device. Press K
% to see a list of all the keyboard commands.
% T. C. O'Haver (toh@umd.edu),  Version 2, October, 2011.  Added : 'Y' key
% cycles through 4 plot modes: linear x,y. linear x, log y; log x, linear
% y; and log-log.  'H' key switchs between frequency and period x scale on
% power spectrum graph in lower window

global samplingtime samplerate f1 f2
global signaltype plotmode
warning off all

% Initial values of parameters
samplingtime=.1;  % Total duration of sampled signal, in sec, millisec, or microsec.
samplerate=8000; % Sample rate in Hz, KHz, or MHz, respectively.
f1=250;
f2=12;
signaltype=1;  % Sine wave initially
plotmode=1;

% Plot the signal and its power spectrum
% set(gcf,'doublebuffer','on'); 
RedrawSpectrum;
figure(1);
% h2=gca;

% Attaches KeyPress test function to the figure.
set(gcf,'KeyPressFcn',@ReadKey)
uicontrol('Style','text')
% end of outer function
% ----------------------------SUBFUNCTIONS--------------------------------

function ReadKey (obj,eventdata)
% Interprets key presses from the Figure window.
% When a key is pressed, interprets key and calls corresponding function.
% Note: If you don't like my key assignments, you can change the numbers
% in the case statements here to re-assign that function to any other key.

global x y samplingtime samplerate f1 f2
global signaltype plotmode XMODE

x=0:(1/samplerate):samplingtime;
key=get(gcf,'CurrentCharacter');
if ischar(key),
  switch double(key),
    case 96
        % '`' key pressed.
        signaltype=1; % Sine wave of frequency F1 (Hz) and phase F2
        RedrawSpectrum;
    case 49
        % 1 key pressed.
        signaltype=13; % Square wave of frequency F1 (Hz) and phase F2;
        RedrawSpectrum;
    case 50
        % 2 key pressed.
        signaltype=3;  % Sawtooth wave of frequency F1(Hz)
        RedrawSpectrum;
    case 51
        %  3 key pressed.
        signaltype=4; % Triangle wave of frequency F1 (Hz) and phase F2
        RedrawSpectrum;
    case 52
        %  4 key pressed.
        signaltype=2; % Sine wave burst of frequency F1 (Hz) and length F2
        RedrawSpectrum;
    case 53
        %  5 key pressed.
        signaltype=5; % 440 Hz carrier amplitude modulated by sine wave of frequency F1 (Hz) and amplitude F2
        RedrawSpectrum;
    case 54
        %  6 key pressed.
        signaltype=6; % 440 Hz carrier frequency modulated by sine wave of frequency F1 (Hz) and amplitude F2
        RedrawSpectrum;
    case 55
        %  7 key pressed.
        signaltype=7; % Sine wave of frequency F1 (Hz) modulated by Gaussian pulse of width F2
        RedrawSpectrum;
    case 56
        %  8 key pressed.
        signaltype=8; % Sine wave of frequency F1 (Hz) with non-linear transfer function F2
        RedrawSpectrum;
    case 57
        %  9 key pressed.
        signaltype=10; % Sine wave sweep from 0 to f1 (Hz)
        RedrawSpectrum;
    case 48
        %  0 key pressed.
        signaltype=9; % Sine wave of frequency F1 (Hz) and amplitude F2 plus random white noise
        RedrawSpectrum;
    case 45
        %  '-' key pressed.
        signaltype=11; % Pink (1/f) noise
        RedrawSpectrum;
    case 61
        %  '=' key pressed.
        signaltype=12; % Sine wave of frequency F1 (Hz) and amplitude F2 plus pink (1/f) noise
        RedrawSpectrum;
    case 8
        % Backspace key pressed.
        signaltype=14; % Used defined signal
        RedrawSpectrum;  
    case 108
        % L key pressed.
        plotmode=plotmode+1;
        if plotmode==5,plotmode=1;end
        RedrawSpectrum;       
    case 104
        % When 'H' is pressed, toggles between frequency and period x=vxis on spectrum 
         if XMODE==0,XMODE=1;else XMODE=0;end
        RedrawSpectrum;     
    case 97
        % When 'a' key is pressed, increases "samplingtime" by 10%
        samplingtime=abs(samplingtime+.1*samplingtime);
        RedrawSpectrum;      
    case 122
        % When 'z' key is pressed, decreases "samplingtime" by 10%
        samplingtime=abs(samplingtime-.1*samplingtime);
        RedrawSpectrum;        
   case 115 % When 's' key is pressed, increases "samplerate" by 10%
         samplerate=samplerate+.1*samplerate;
         RedrawSpectrum;    
   case 120 % When 'x' key is pressed, decreases "samplerate" by 10%
         samplerate=samplerate-.1*samplerate;
         RedrawSpectrum;    
   case 100
        % When 'd' key is pressed, increases "f1" by 10%
        f1=abs(f1+.1*f1);
        RedrawSpectrum;   
    case 99
        % When 'c' key is pressed, decreases "f1" by 10%
        f1=abs(f1-.1*f1);
        RedrawSpectrum;   
    case 102
        % When 'f' key is pressed, increases "f2" by 10%
         f2=abs(f2+.1*f2);
         RedrawSpectrum;    
    case 118
        % When 'v' key is pressed, decreases "f2" by 10%
         f2=abs(f2-.1*f2);
         RedrawSpectrum; 
    case 112
        % When 'p' key is pressed
        soundsc(y./2,samplerate)
    case 13
        % When 'Enter' key is pressed
        soundsc(y./2,samplerate)
    case 107
        % When 'k' key is pressed, prints out table of keyboard commands
        disp('                                                 ')
        disp('KEYBOARD CONTROLS:')
        disp(' Adjust signal duration 10% up/down.........A,Z')
        disp(' Adjust sampling rate 10% up/down...........S,X')
        disp(' Adjust first variable 10% up/down......... D,C')
        disp(' Adjust second variable 10% up/down........ F,V')
        disp(' Cycle through Linear/Log plot modes..........L')
        disp(' Switch X-axis scale of power spectrum........H')
        disp(' Print keyboard commands......................K')
        disp(' Play signal as sound................Enter or P')
        disp('                                                 ')
        disp('SIGNAL TYPES, selected by  ` 1 2 3 ... - = Backspace')
        disp(' Sine wave of frequency F1 (Hz) ') 
        disp('     and phase F2.............................. `')
        disp(' Square wave of frequency F1 (Hz) ') 
        disp('     and phase F2...............................1')
        disp(' Sawtooth wave of ') 
        disp('      frequency F1(Hz)..........................2')
        disp(' Triangle wave of frequency F1 (Hz) ') 
        disp('     and phase F2...............................3')
        disp(' Sine wave burst of frequency F1 (Hz) ') 
        disp('     and length F2..............................4')
        disp(' 440 Hz carrier amplitude modulated by sine') 
        disp('  wave of frequency F1 (Hz) and amplitude F2....5')
        disp(' 440 Hz carrier frequency modulated by sine') 
        disp('  wave of frequency F1 (Hz) and amplitude F2....6')
        disp(' Sine wave of frequency F1 (Hz) modulated') 
        disp('     with Gaussian of width F2 sec..............7')
        disp(' Sine wave of frequency F1 (Hz) with non-') 
        disp('     linear transfer function F2................8')
        disp(' Sine wave sweep from 0 to f1 (Hz)..............9')
        disp(' Sine wave of frequency F1 (Hz) and amplitude') 
        disp('      F2 plus random white noise................0')
        disp(' Pink (1/f) noise...............................-')
        disp(' Sine wave of frequency F1 (Hz) and amplitude') 
        disp('      F2 plus pink noise........................=')
        disp(' Used-defined signal waveform...........Backspace')
      otherwise  
       UnassignedKey=double(key)
       disp('Press k to print out list of keyboard commands')
   end % switch
end % if
% ----------------------------------------------------------------------

function RedrawSpectrum
% Re-computes amd plots the signal and the power 
% spectrum each time one of the keys is pressed.

global x y samplingtime samplerate f1 f2
global signaltype signalstring plotmode XMODE

x=0:(1/samplerate):samplingtime;

% Computes signal for the selected signal type.
switch signaltype
    case 1
        signalstring=sprintf('Sine wave, frequency %.3g Hz (D,C to adjust), phase %.2g (F,V to adjust)',f1,f2);
        y=sin(f1*2*pi.*(x+(f2-1)/10000));
    case 2 ;
        signalstring=sprintf('Sine wave burst, frequency %.3g Hz (D,C to adjust), length %.2g sec (F,V to adjust)',f1,2*f2/samplerate);
        y=sin(f1*2*pi.*x);
        m=round(length(y)/2);
        if f2<1,f2=1;end
        y(1:round(m-f2))=0;
        y(round(m+f2):length(y))=0;
    case 3  
        signalstring=sprintf('Sawtooth wave of frequency %.3g Hz (D,C to adjust)',f1);
        y=rem(x,max(x)/((f1)*samplingtime));y=y-mean(y);y=y/max(y);
    case 4 
        signalstring=sprintf('Triangle wave, frequency %.3g Hz (D,C to adjust), phase %.2g (F,V to adjust)',f1,f2);
        y=asin(sin(f1*2*pi.*(x+(f2-1)/10000)))/1.568;
    case 5 
        signalstring=sprintf('440 Hz carrier AM by %.3g Hz (D,C to adjust), amplitude %.2g (F,V to adjust)',f1,f2);
        y=(1+f2/100.*sin(2*f1*pi.*x)).*sin(880*pi.*x); 
    case 6
        signalstring=sprintf('440 Hz carrier FM by sine wave %.3g Hz (D,C to adjust), amplitude %.2g (F,V to adjust)',f1,f2);
        y=sin(880*pi*x+f2/50.*sin(2*pi*f1*x)); 
    case 7
        signalstring=sprintf('%.3g Hz sine wave (D,C to adjust) modulated by a Gaussian of width %.2g sec (F,V to adjust)',f1,2*f2*max(x)/200);
        y=sin(2*f1*pi.*x).*gaussian(x,max(x)/2,f2*max(x)/100);
    case 8
        signalstring=sprintf('Sine wave, frequency %.3g Hz (D,C to adjust), with non-linear transfer function %.2g (F,V to adjust)',f1,f2);
        y=(sin(f1*2*pi.*(x+(f2-1)/10000)));
        y=real(y.^(1+f2/100));
    case 9
        signalstring=sprintf('White noise plus sine wave, frequency %.3g Hz (D,C to adjust), amplitude %.2g (F,V to adjust)',f1,f2);
        y=f2/10*sin(f1*2*pi.*x)+randn(size(x));y=y-mean(y);y=y/max(y);  
    case 10   
        signalstring=sprintf('Sine wave sweep from 0 to %.3g Hz (D,C to adjust)',f1);
        y=sin(x/max(x)*f1*2*pi.*(x));
    case 11
        signalstring='Pink (1/f) noise';
        y=pinknoise(length(x)); y=y-mean(y);y=y/max(y);     
    case 12
        signalstring=sprintf('Pink (1/f) noise p1us sine wave, frequency %.3g Hz (D,C to adjust), amplitude %.2g (F,V to adjust)',f1,f2/10);
        y=f2/10*sin(f1*2*pi.*x)+pinknoise(length(x));y=y-mean(y);
        % y=y/max(y);  
    case 13
        signalstring=sprintf('Square wave, frequency %.3g Hz (D,C to adjust), phase %.2g (F,V to adjust)',f1,f2);
        y=sign(sin((f1+f2/100)*2*pi.*(x+(f2-1)/1000))); 
    case 14
        signalstring='User-defined signal (change in line 254)';
        y=sin(f1*2*pi.*(x+(f2-1)/10000)); %  Replace with your own function (x=time, 
end % switch signaltype

% Plot the signal in the upper half of the window.
subplot(2,1,1)
plot(x,y)
title( num2str(signalstring) )
xlabel([ 'Signal duration = ' num2str(0.001*round(1000*samplingtime)) 'sec (A,Z to adjust)   Sampling rate = ' num2str(round(samplerate)) 'Hz (S,X to adjust)' ])
axis([min(x) max(x) min(y) max(y)])

% Plot the power spectrum  in the lower half of the window.
subplot(2,1,2)
fy=fft(y);
py=sqrt(fy .* conj(fy)); % Compute power spectrum
plotrange=1:length(fy)/2;

if XMODE,
   f=range(x)./(plotrange-1);
else
   f=(plotrange-1)./range(x);
end

switch plotmode,
  case 1
    plot(f,real(py(plotrange)),'r')
    ylabel('Linear y')
  case 2
    semilogx(f,real(py(plotrange)),'r')
    ylabel('Linear y')
  case 3
    semilogy(f,real(py(plotrange)),'r')
    ylabel('Log y')
  case 4
    loglog(f,real(py(plotrange)),'r')
    ylabel('Log y')
    otherwise,
end
title('Power spectrum of the above signal. Press K  to display key commands.')
if XMODE,
    xlabel('x=Time (e.g. seconds)  Press H to change to Frequency')
else
    xlabel('x=Frequency (e.g. cycles/second}  Press H to change to Time')
end
axis([min(f) max(f) 0 max(py)])
% ---------------------------------------------------------------------- 

function ry=pinknoise(n)
% Random noise with pink (1/f) power spectrum with mean zero 
% and unit standard deviation. n is number of points.
% Tom O'Haver, 2008
x=1:n;
y=randn(size(x));  % Random normally-distributed white noise
% Fourier filter 
fy=fft(y); % Compute Fourier transform of signal y
% Compute filter shape
lft1=1:(length(fy)/2)+1;
lft2=(length(fy)/2):length(fy);
ffilter1=ones(size(lft1))./(sqrt(lft1));
ffilter2=ones(size(lft2))./(sqrt(lft2));
ffilter=[ffilter1,ffilter2];
if length(fy)>length(ffilter), ffilter=[ffilter ffilter(1)];end
ffy=fy.*ffilter(1:length(fy));  % Multiply filter by Fourier transform of signal
ry=real(ifft(ffy)); % Inverse transform to recover filtered signal 'ry'
ry=((ry-mean(ry))./std(ry)); % Normalize to unit standard deviation
% ----------------------------------------------------------------------    
function g = gaussian(x,pos,wid)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Example: gaussian([1 2 3],1,2) gives result [0.5000    1.0000    0.5000]
g = exp(-((x-pos)./(0.6005612.*wid)) .^2);

function r = range(arr)
r = max(arr) - min(arr);
