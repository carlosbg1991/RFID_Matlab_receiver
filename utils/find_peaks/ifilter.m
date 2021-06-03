function ry=ifilter(ix,iy,icenter,iwidth,ishape,imode,ifilt)
% ifilter(x,y) or ifilter(y) or ifilter([x y]) or
% ry=ifilter(x,y,center,width,shape,plotmode,filtermode)
% Keyboard-operated interactive Fourier filter function for
% time-series signal (x,y), with keyboard controls that 
% allow you to adjust the filter parameters continuously 
% while observing the effect on your signal dynamically. Optinal
% input arguments set the intital values of center frequency, 
% filter width, shape, plotmode (1=linear; 2=semilog frequency; 
% 3=semilog amplitude; 4=log-log) and filtermode ('Band-pass',
% 'Lowpass', 'Highpass', 'Bandpass', and 'Band-reject (notch)')
% 'X' key changes x-axis scale on spectrum plot between frequency
% and time (period); zeroth harmonic (DC level) skipped in spectrum plot.
% Returns the filtered signal. Press K to list keyboard commands.
% T. C. O'Haver (toh@umd.edu), Version 4.3, June, 2016. 
%
% Example 1:
% Periodic waveform with 3 harmonic components
% x=0:1000; y=sin(x/10)+sin(x/5)+sin(x);ifilter(x,y);
%
% Example 2: uses optional input arguments to set initial values:
% x=0:(1/8000):.3;
% y=(1+12/100.*sin(2*47*pi.*x)).*sin(880*pi.*x)+(1+12/100.*sin(2*20*pi.*x)).*sin(2000*pi.*x);
% ry=ifilter(x,y,440,31,18,3,'Band-pass');
%
% Example 3: Picking one frequency out of a noisy sine wave.
% x=[0:.01:2*pi]';
% y=sin(20*x)+3.*randn(size(x));
% ifilter(x,y,3.1,0.85924,15,1,'Band-pass');
%
% Example 4: Square wave with band-pass vs Comb pass filter
% t = 0:.0001:.0625;
% y=square(2*pi*64*t);
% ifilter(t,y,64,32,12,1,'Band-pass');
% ifilter(t,y,48,32,2,1,'Comb pass');
%
% KEYBOARD CONTROLS when figure window is topmost:
% Adjust center frequency.......Coarse: < and >
%                               Fine: left and right cursor arrows
% Adjust filter width...........Coarse: / and "  
%                               Fine: up and down cursor arrows
% Filter shape..................A,Z (A more rectangular, Z more Gaussian)
% Filter mode...................B=bandpass; N or R=notch (band reject)
%                               H=High-pass; L=Low-pass
% Select plot mode..............B=bandpass; N or R=notch (band reject);H=High-pass;
%                               L=Low-pass; C=Comb pass; V=Comb notch
% Print keyboard commands.......K  Pints this list
% Print filter parameters.......Q  Prints input arguments: center,width,shape,plotmode,filtermode
% Print current settings........T  Prints list of current settings
% Switch SPECTRUM X-axis scale..X switch between frequency and period x scale on POWER SPECTRA
% Switch OUTPUT Y-axis scale....Y switch between fixed or variable y scale on output plot
% Play output as sound..........P or Enter
% Save output as .mat file......S
%  
global x y CENTER WIDTH SHAPE PLOTMODE FILTERMODE YMODE XMODE
format short g
format compact
switch nargin
    % 'nargin' is the number of arguments
    case 1
        datasize=size(ix);
        if isvector(ix),
            x=1:length(ix); % Use this only to create an x vector if needed
            y=ix;
        else
            if datasize(1)<datasize(2),ix=ix';end
            x=ix(:,1);
            y=ix(:,2);
        end
        % Adjust x and y vector shape to 1 x n (rather than n x 1)
        x=reshape(x,1,length(x));
        y=reshape(y,1,length(y));
        % If necessary, flip the data vectors so that x increases
        if x(1)>x(length(x)),
            disp('x-axis flipped.')
            x=fliplr(x);
            y=fliplr(y);
        end
        CENTER=1;
        WIDTH=10;
        SHAPE=2;
        FILTERMODE='Band-pass'; % Starts in band-pass mode
        PLOTMODE=1;  % Starts in plotmode=1 for linear plot
    case 2
        x=ix;
        y=iy;
        % Set initial values of filter parameters.
        CENTER=1;
        WIDTH=10;
        SHAPE=2;
        FILTERMODE='Band-pass'; % Starts in band-pass mode
        PLOTMODE=1;  % Starts in plotmode=1 for linear plot
    case 7
        x=ix;
        y=iy;
        CENTER=icenter.*range(x);
        WIDTH=iwidth.*range(x);
        SHAPE=ishape;
        FILTERMODE=ifilt;
        PLOTMODE=imode;
    otherwise
        disp('Invalid number of arguments')
end % switch nargin

ry=ones(size(x));
YMODE=0;XMODE=0;
% Adjust x and y vector shape to 1 x n (rather than n x 1)
x=reshape(x,1,length(x));
y=reshape(y,1,length(y));

% Plot the signal and its power spectrum
figure(1)
ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);

% Attaches KeyPress test function to the figure.
set(gcf,'KeyPressFcn',@ReadKey)
uicontrol('Style','text')
% end of outer function
% ----------------------------SUBFUNCTIONS--------------------------------
function ReadKey(obj,eventdata)
% Interprets key presses from the Figure window.
% When a key is pressed, interprets key and calls corresponding function.
% Note: If you don't like my key assignments, you can change the numbers
% in the case statements here to re-assign that function to any other key.
global x y CENTER WIDTH SHAPE PLOTMODE FILTERMODE ry YMODE XMODE
key=get(gcf,'CurrentCharacter');
if ischar(key),
  switch double(key),
    case 28
        % Pans "CENTER" one point down when left arrow pressed.
          CENTER=CENTER-CENTER/100;
          if CENTER<1,CENTER=1;end
          ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 29
        % Pans "CENTER" one point up when left arrow pressed.
          CENTER=CENTER+CENTER/100;
          ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 44
        % Pans "CENTER" 10 points down when < key pressed.
          CENTER=CENTER-2*WIDTH;
          if CENTER<1,CENTER=1;end
          ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 46
        % Pans "CENTER" WIDTH points up when > key pressed.
          CENTER=CENTER+2*WIDTH;
          ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 30
        % Zooms filter width one point up when up arrow pressed.
          WIDTH=WIDTH+WIDTH/100;
          ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 31
        % Zooms filter width one point down when down arrow pressed.
          WIDTH=WIDTH-WIDTH/100;
          if WIDTH<1,WIDTH=1;end
          ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 39
        % Zooms filter width up when / pressed.
          WIDTH=WIDTH+WIDTH/10;
          ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 47
        % Zooms filter width down when ' pressed.
          WIDTH=WIDTH-WIDTH/10;
          if WIDTH<1,WIDTH=1;end
          ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 122
        % When 'z' key is pressed, shape is made more rectangular
          SHAPE=SHAPE-1;
          if SHAPE<0,SHAPE=0;end
          ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 97
        % When 'a' key is pressed, shape is made more Gaussian
          SHAPE=SHAPE+1;
          ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 49 % When '1' key is pressed, set plot mode to linear/linear
         PLOTMODE=1;
         ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 50 % When '2' key is pressed, set plot mode to log x, linear y
         PLOTMODE=2; 
         ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 51
        % When '3' key is pressed, set plot mode to linear x, log y
        PLOTMODE=3;
        ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 52
        % When '4' key is pressed,  set plot mode to log y, log x
        PLOTMODE=4;
        ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 98
        % When 'b' key is pressed, set to bandpass mode
         FILTERMODE='Band-pass';
         ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case {110, 114}
        % When 'n' or 'r' key is pressed, set to band reject (FILTERMODE) node
         FILTERMODE='Band-reject (notch)';
         ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 104
        % When 'H' key is pressed, makes high-pass filter
        FILTERMODE='High-pass';
        ry=RedrawFourierFilter(x,y,length(y)/2,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 108
          % When 'L' is pressed, makes low=pass filter
          FILTERMODE='Low-pass';
          ry=RedrawFourierFilter(x,y,0,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 99
          % When 'C' is pressed, makes comb filter
          FILTERMODE='Comb pass';
          ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 118
          % When 'V' is pressed, makes comb filter
          FILTERMODE='Comb notch';
          ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
    case 13
        % When 'Enter' is pressed, plays filtered signal as sound
        % without scaling
        sound(ry,8000);
    case {112, 32}
        % When 'p' key is pressed, plays filtered signal as sound, scaled 
        % so that the sound is played as loud as possible without clipping.
          sound(ry./max(ry),44000)  
    case 115
        % When 's' key is pressed, filtered signal is saved as
        % FilteredOutput.mat
        save FilteredOutput ry
    case 107
        % When 'k' key is pressed, prints out table of keyboard commands
        disp('')   
        disp('KEYBOARD CONTROLS when figure window is topmost:')
        disp('Adjust center frequency.......Coarse: < and >')   
        disp('                              Fine: left and right cursor arrows')
        disp('Adjust filter width...........Coarse: / and "  ') 
        disp('                              Fine: up and down cursor arrows')
        disp('Filter shape..................A,Z (A more rectangular, Z more Gaussian)')
        disp('Filter mode...................B=bandpass; N or R=notch (band reject); H=High-pass;')
        disp('                              L=Low-pass; C=Comb pass; V=Comb notch (reject)')        
        disp('Select plot mode..............1=linear; 2=semilog frequency' )
        disp('                              3=semilog amplitude; 4=log-log')
        disp('Print keyboard commands.......K  Pints this list')
        disp('Print filter parameters.......Q  Prints input arguments: center,width,shape,plotmode,filtermode')
        disp('Print current settings........T  Prints list of current settings')
        disp('Switch SPECTRUM X-axis scale..X switch between frequency and period x scale on POWER SPECTRA')
        disp('Switch OUTPUT Y-axis scale....Y switch between fixed or variable y scale on output plot')
        disp('Play output as sound..........P or Enter')
        disp('Save output as .mat file......S')
        disp(' ')     
      case 113
        % When 'Q' is pressed, prints fourier filter parameters on a single line     
        disp([ 'ifilter(x,y,' num2str(CENTER./range(x)) ',' num2str(WIDTH./range(x))  ',' num2str(SHAPE)  ',' num2str(PLOTMODE) ',''' num2str(FILTERMODE) ''');'])
      case 119
        % When 'W' is pressed, prints fourier filter parameters on a single line     
        disp([ 'ifilter(x,y,' num2str(CENTER./range(x)) ',' num2str(WIDTH./range(x))  ',' num2str(SHAPE)  ',' num2str(PLOTMODE) ',''' num2str(FILTERMODE) ''');']) 
      case 116
        % When 'T' is pressed, prints list of current settings
       disp('------ iFilter Settings --------' )
       disp(['Number of points in signal: ' num2str(length(y)) ] )
       disp(['Duration of signal (x range): ' num2str(range(x)) ] )
       disp(['Interval between x-axis values: ' num2str(x(2)-x(1)) ] )
       disp(['Center harmonic: ' num2str(CENTER) ] )
       disp(['Center frequency: ' num2str(CENTER./range(x)) ] ) 
       disp(['Center period: ' num2str(range(x)./CENTER) ] )       
       disp(['Width in harmonics: ' num2str(WIDTH) ] )
       disp(['Width in frequency: ' num2str(WIDTH./range(x)) ] )  
       disp(['Filter mode: ' FILTERMODE ] ) 
       disp('To call iFilter with these parameters, type:') 
       disp(['ifilter(x,y,' num2str(CENTER./range(x)) ',' num2str(WIDTH./range(x))  ',' num2str(SHAPE)  ',' num2str(PLOTMODE) ',''' num2str(FILTERMODE) ''');']) 
      case 121
         % When 'Y' is pressed, toggles between fixed or variable y scale on OUTPUT plot 
         if YMODE==0,YMODE=1;else YMODE=0;end
         ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
      case 120
         % When 'X' is pressed, toggles between frequency and period x=vxis on spectrum 
         if XMODE==0,XMODE=1;else XMODE=0;end
         ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
      otherwise  % if key pressed is unnasigned
       UnassignedKey=double(key) % print unassigned key code
       disp('Press K to print out list of keyboard commands')
   end % switch
end % if

function ry=RedrawFourierFilter(xvector,yvector,centerfrequency,filterwidth,filtershape,mode,FILTERMODE,YMODE,XMODE)
% Separate graph windows for the original and filtered signals.
% Computes and plots fourier filter for signal yvector.  
% Centerfrequency and filterwidth are the center frequency and
% width of the pass band, in harmonics, 'filtershape' determines
% the sharpness of the cut-off. Plot modes: mode=1 linear x and y; 
% mode=2 log x linear y; mode=3 linear x; log y; mode=4 log y log x
fy=fft(yvector);
lft1=[1:(length(fy)/2)];
lft2=[(length(fy)/2+1):length(fy)];
% Compute filter shape.
if strcmp(FILTERMODE,'Band-pass'),
       ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
       ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
       ffilter=[ffilter1,ffilter2];
end
if strcmp(FILTERMODE,'High-pass')
       centerfrequency=length(xvector)/2;
       ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
       ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
       ffilter=[ffilter1,ffilter2];
end
if strcmp(FILTERMODE,'Low-pass')
       centerfrequency=0;
       ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
       ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
       ffilter=[ffilter1,ffilter2];
end
if strcmp(FILTERMODE,'Band-reject (notch)')
       ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
       ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
       ffilter=1-[ffilter1,ffilter2]; 
end
if strcmp(FILTERMODE,'Comb pass')
    n=2;
    ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
    ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
    while n<50,
       ffilter1=ffilter1+shape(lft1,n*(centerfrequency+1),filterwidth,filtershape);
       ffilter2=ffilter2+shape(lft2,length(fy)-n*(centerfrequency+1),filterwidth,filtershape);
       n=n+1;
    end
       ffilter=[ffilter1,ffilter2];
end
if strcmp(FILTERMODE,'Comb notch')
    n=2;
    ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
    ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
    while n<30,
       ffilter1=ffilter1+shape(lft1,n*(centerfrequency+1),filterwidth,filtershape);
       ffilter2=ffilter2+shape(lft2,length(fy)-n*(centerfrequency+1),filterwidth,filtershape);
       n=n+1;
    end
       ffilter=1-[ffilter1,ffilter2];
end
if length(fy)>length(ffilter), ffilter=[ffilter ffilter(1)];end
ffy=fy.*ffilter;  % Multiply filter by Fourier Transform of signal
ry=real(ifft(ffy));

subplot(3,1,1)  % Plot original signal in top plot
plot(xvector,yvector,'b');
title('iFilter 4.3:       Press K to display keyboard commands.')
xlabel('INPUT: Original Signal   x=time')
axis([xvector(1) xvector(length(xvector)) min(yvector) max(yvector)]);  

subplot(3,1,3)  % Plot filtered signal in lower plot
plot(xvector,ry,'r');
if YMODE,
   axis([xvector(1) xvector(length(xvector)) min(ry) max(ry)]);  
else
  axis([xvector(1) xvector(length(xvector)) min(yvector) max(yvector)]);    
end
title([ FILTERMODE ' mode:  Freq= ' num2str(centerfrequency./range(xvector)) '    Period= ' num2str(range(xvector)./centerfrequency) '   Width= ' num2str(filterwidth./range(xvector)) '     Shape=  ' num2str(filtershape)])
xlabel('OUTPUT: Filtered signal    x=time')
if YMODE==0,
    ylabel('Fixed Y scale')
else
    ylabel('Variable Y scale')
end
subplot(3,1,2)    % Plot power spectrum and filter in middle plot
py=fy .* conj(fy); % Compute power spectrum
plotrange=2:length(fy)/2;
if XMODE,
   f=range(xvector)./(plotrange-1);
else
   f=((plotrange-1)./range(xvector));
end
switch mode,
  case 1
    plot(f,real(py(plotrange)),f,max(real(py(plotrange))).*ffilter(plotrange),'r')
    ylabel('Linear y')
  case 2
    semilogx(f,real(py(plotrange)),f,max(real(py(plotrange))).*ffilter(plotrange),'r')
    ylabel('Linear y')
  case 3
    semilogy(f,real(py(plotrange)),f,max(real(py(plotrange))).*ffilter(plotrange),'r')
    ylabel('Log y')
  case 4
    loglog(f,real(py(plotrange)),f,max(real(py(plotrange))).*ffilter(plotrange),'r')
    ylabel('Log y')
    otherwise,
end
title('POWER SPECTRA:  BLUE = Input signal    RED = Filter')
if XMODE,
    xlabel('x=Time (e.g. seconds)')
else
    xlabel('x=Frequency (e.g. cycles/second)')
end
% Expand frequency axis if filter is in the lower half of frequency axis
if centerfrequency+filterwidth/2<length(fy)/4,
    axis([0 max(f)/2 min(py) 1.1*max(real(py(plotrange)))]) 
else
    axis([0 max(f) min(py) 1.1*max(real(py(plotrange)))])
end

function g = shape(x,pos,wid,n)
%  shape(x,pos,wid,n) = peak centered on x=pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  Shape is Lorentzian (1/x^2) when n=0, Gaussian (exp(-x^2))
%  when n=1, and becomes more rectangular as n increases.
%  Example: shape([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
if n==0
    g=ones(size(x))./(1+((x-pos)./(0.5.*wid)).^2);
else
    g = exp(-((x-pos)./(0.6.*wid)) .^(2*round(n)));
end

function r = range(arr)
r = max(arr) - min(arr);