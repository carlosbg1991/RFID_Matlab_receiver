function idemo
% Self-contained demonstration function for comparing the ipeak and
% peakfit functions applied to a test signal consisting of several
% narrow peaks and one broad peak in the center. 
%  T. C. O'Haver, August 2012

increment=1;
x=300:increment:2000;

% For each simulated peak, compute the amplitude, position, and width
pos=[400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900];   % Positions of the peaks (Change if desired)
amp=round(1+5.*rand(1,length(pos)));  % Amplitudes of the peaks  (Change if desired)
wid=[5 5 5 5 5 5 5 200 5 5 5 5 5 5 5 5];
amp(16)=10;
Noise=0.05; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(k,:) = [k pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
demodata=[x' y']; % Assembles x and y vectors into data matrix

disp('-----------------------------------------------------------------')
disp('Demonstration function for comparing the ipeak.m and peakfit.m ')
disp('functions applied to a test signal consisting of several')
disp('narrow spikes and one broad peak in the center.')
disp(' ')
disp('Using iPeak:')
disp('iPeak can be adjusted to detect some peaks and ignore others')
disp('by means of the input arguments or keystroke controls.')
disp(' ')
disp('>> iPeakResults=ipeak(demodata,0,AmpT,SlopeT,SmoothW,FitW); ')
disp(' ')
disp('AmpT - Discriminates on the basis of peak height. Any peaks with')
disp(' height less than this value are ignored. Normally this is set to')
disp(' lower than the smallest peak but higher than the noise.')
disp(' ')
disp('SlopeT - Discriminates on the basis of peak width. Larger ')
disp(' values of this parameter will neglect broad features of the')
disp(' signal. A reasonable initial value for Gaussian peaks is 0.7*W^-2,')
disp(' where W is the number of data points in the half-width of the peak.')
disp(' ')
disp('SmoothW - Width of the smooth function that is applied to data before')
disp(' the slope is measure d. Larger values of SmoothW will neglect small,')
disp(' sharp features. A reasonable value is typically about equal to 1/2 of')
disp(' the number of data points in the half-width of the peaks.')
disp(' ')
disp('FitW - The number of points around the "top part" of the (unsmoothed)')
disp(' peak that are taken to estimate the peak heights, positions, and widths.')
disp(' A reasonable value is typically about equal to 1/2 of the number of')
disp(' data points in the half-width of the peaks. The minimum value is 3.')
disp(' ')
fprintf(2,'Press any key to continue.....or press Ctrl-C to abort.\n')
pause
disp('-----------------------------------------------------------------')
disp('For example, to detect the broader peaks and ignore the very narrow')
disp('peaks and spikes, use a small value of SlopeT and a large value of')
disp('SmoothW and FitW: ')
disp(' ')
disp('>> AmpT=0.2')
disp('>> SlopeT=0.00001')
disp('>> SmoothW=121')
disp('>> FitW=121')
disp('>> iPeakResults=ipeak(demodata,0,AmpT,SlopeT,SmoothW,FitW)')
disp(' ')
% Now call iPeak, with specified values of the peak detection parameters:
% AmpT, SlopeT, SmoothW, and FitW;
tic;
iPeakResults=ipeak(demodata,0,.2,0.00001,121,121,1100,600);
iPeakTime=toc;
NumPeaks=max(iPeakResults(:,1));
disp(['Number of peaks detected = ' num2str(NumPeaks)]) 
disp(['Elapsed time is ' num2str(iPeakTime) ' seconds'])
disp('         Peak #    Position      Height       Width        Area')
iPeakResults
drawnow
disp(' ')
disp('In this case, iPeak detects and measures only the broad peak at')
disp('x=1100 (whose real width was 200 units before noise was added).')
disp(' ')
fprintf(2,'Press any key to continue.....or press Ctrl-C to abort.\n')
pause
disp('-----------------------------------------------------------------')
disp('On the other hand, to detect the narrow peaks and ignore the ')
disp('broad peaks, use a larger value of SlopeT and a small value of')
disp('SmoothW and FitW: ')
disp(' ')
disp('>> AmpT=0.2')
disp('>> SlopeT=0.01')
disp('>> SmoothW=7')
disp('>> FitW=6')
disp('>> iPeakResults=ipeak(demodata,0,AmpT,SlopeT,SmoothW,FitW)')
tic;
iPeakResults=ipeak(demodata,0,.2,0.01,7,6,400,50);
iPeakTime=toc;
NumPeaks=max(iPeakResults(:,1));
disp(' ')
disp(['Number of peaks detected = ' num2str(NumPeaks)]) 
disp(['Elapsed time is ' num2str(iPeakTime) ' seconds'])
disp('         Peak #    Position      Height       Width        Area')
iPeakResults 
disp(' ')
disp('In this case, iPeak detects and measures only the narrow peaks')
disp('and skips the broad peak at 1100.') 
disp('  You can easily estimate the accuracy of the measurement of peak ')
disp('position, height, and width because all the peaks have true peak')
disp('positions that are exactly on the 100s (400, 500, 600, etc),')
disp('integer peak heights (1,2,3...6), and peak widths of either 200 or 5.')
disp('  The accuracy of peak position is the best. The accuracy of peak')
disp('and width depend on the situation: isolated peaks of Gaussian shape')
disp('are measured accurately, but peaks on a background, overlapping peaks,')
disp('and non-Gaussian peaks require special care. For the best results')
disp('in those cases, it''s best to use iterative peak fitting, which')
disp('can be done within iPeak or separately using the peakfit.m ')
disp('or ipf.m functions. For example.....')
disp(' ')
fprintf(2,'Press any key to continue.....or press Ctrl-C to abort.\n')
pause
disp('-----------------------------------------------------------------')

disp('The narrow peak at x=1000 is a good example.  Its true height')
disp(['is ' num2str(ActualPeaks(7,3)) ' and its width is ' num2str(ActualPeaks(7,4)) ', but iPeak did not measure those' ] )
disp('values accurately because that peak sits on the side of the')
disp('broad peak at x=1100, which acts as a background under the peak')
disp('at 1000, making both the peak height and with too high:') 
disp('  ') 
disp('         Peak #    Position      Height       Width        Area')
iPeak_Results=iPeakResults(7,:)
disp('  ')
disp('There are two ways to handle this: (a) fit the two peaks together')
disp('using peakfit.m, or (b) use the "autozero" mode to correct ')
disp('the peak at 1000 for the background near that peak.')
disp('  ')
disp('(a) Using peakfit.m and a two-peak fit: ')
disp('Here we treat the background as another peak to be fit.')
disp(' >> FitResults=peakfit(demodata,1055,225,2,1,0,10,0,0)')
disp('         Peak #    Position      Height       Width        Area')
FitResults=peakfit(demodata,1055,225,2,1,0,10,0,0) 
disp('  ')
disp('The peak at 1000 (Peak 1 in this table) is measured much more ')
disp('accurately here than iPeak did. The broader peak at 1100 is also')
disp('measured, but you may not be interested in that peak.')
disp('  ')
fprintf(2,'Press any key to continue.....or press Ctrl-C to abort.\n')
pause
disp('-----------------------------------------------------------------')
disp('(b) Using Autozero and a one-peak fit: ')
disp('Both iPeak and peakfit (as well as peakfit''s interactive version,')
disp('ipf.m) all have an Autozero mode that subtracts a linear background')
disp('from under a peak, assuming that the background close to the peak')
disp('is locally linear (or nearly so).')
disp('  ')
disp('Using iPeak with Autozero turned on (9th input argument=1):')
disp('>> iPeakResults=ipeak(demodata,0,.2,0.01,7,6,1000,46,1);  ')
disp('         Peak #    Position      Height       Width        Area')
iPeakResults=ipeak(demodata,0,.2,0.01,7,6,1000,46,1);
iPeakResults(7,:) 
disp('The results for the peak at 1000 (#7) are much better with Autozero on.')
disp('  ')
fprintf(2,'Press any key to continue.....or press Ctrl-C to abort.\n')
pause
disp('  ')
disp('Using peakfit with Autozero turned on (9th input argument=1):')
disp('>> FitResults=peakfit(demodata,1000,46,1,1,0,10,0,1)')
disp('         Peak #    Position      Height       Width        Area')
FitResults=peakfit(demodata,1000,46,1,1,0,10,0,1) 
disp('  ')
disp('The accuracy of measurement is comparable to the previous 2-peak fit.')
disp('  ')
disp('iPeak can also use both these methods: (a) by using the N key to')
disp('perform a peakfit on the peaks numbered in the upper window, and')
disp('(b) by using the T key to toogle on the Autozero mode. See  ')
disp('http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm')
disp('for details.  ')
disp('  ')
disp('End of demo')
disp('If you run this again, you''ll get slightly different results')
disp('because random noise is added to the signal each time. You may')
disp('change the Noise level in line 15 if you wish.')
% -------------------------------------------------------------------

function P=ipeak(DataMatrix,PeakD,AmpT,SlopeT,SmoothW,FitW,xcenter,xrange,MaxError,positions,names)
global X Y xo dx SlopeThreshold AmpThreshold SmoothWidth FitWidth AUTOZERO valleymode
global PeakLabels PeakID Names Positions maxerror logplot plotcolor showpeak

format short g
format compact
warning off all
switch nargin
    % 'nargin' is the number of arguments
    case 1  % One argument only
      % Assumne that the argument must be a matrix of data.
      % If DataMatrix is in the wrong transposition, fix it.
      datasize=size(DataMatrix);
      if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
      datasize=size(DataMatrix);
      if datasize(2)==1, %  Must be ipeak(Y-vector)
         X=[1:length(DataMatrix)]'; % Create an independent variable vector
         Y=DataMatrix;
      else
         % Must be ipeak(DataMatrix)
         X=DataMatrix(:,1); % Split matrix argument 
         Y=DataMatrix(:,2);
      end
      % Calculate default values of peak detection parameters
      PeakDensity=20;   
      % Estimate approx number of points in a peak half-width
      WidthPoints=length(Y)/PeakDensity;  
      SlopeThreshold=WidthPoints^-2;  
      AmpThreshold=abs(min(Y)+0.1*(max(Y)-min(Y))); 
      SmoothWidth=round(WidthPoints/3);  
      FitWidth=round(WidthPoints/3);
      if FitWidth>100,FitWidth=100;end   % Keep FitWidth below 100
      if SmoothWidth>100,SmoothWidth=100;end   % Keep SmoothWidth below 100
      xo=length(Y)/2; % Initial Pan setting
      dx=length(Y)/4; % Initial Zoom setting
      AUTOZERO=0;
    case 2
      % Two arguments; might be separate x and y data vectors, 
      % data matrix + number, or y vector + number (peak density estimate)
      if isscalar(PeakD)    
         % Must be one data matrix and a peak density estimate.
         % If DataMatrix is in the wrong transposition, fix it.
          datasize=size(DataMatrix);
         if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
         datasize=size(DataMatrix);
         if datasize(2)==1, %  Must be ipeak(Y-vector)
            X=[1:length(DataMatrix)]'; % Create an independent variable vector
            Y=DataMatrix;
         else
            % Must be ipeak(DataMatrix)
            X=DataMatrix(:,1); % Split matrix argument 
            Y=DataMatrix(:,2);
         end
         % Calculate values of peak detection parameters
         % arguments based on the peak density, PeakD
         PeakDensity=PeakD;    
         % Estimate approx number of points in a peak half-width
         WidthPoints=length(Y)/PeakDensity;  
         SlopeThreshold=WidthPoints^-2;  
         AmpThreshold=abs(min(Y)+0.1*(max(Y)-min(Y)));  
         SmoothWidth=round(WidthPoints/3);  
         FitWidth=round(WidthPoints/3);
         if FitWidth>100,FitWidth=100;end   % Keep FitWidth below 100
         if SmoothWidth>100,SmoothWidth=100;end   % Keep SmoothWidth below      xo=length(Y)/2; % Initial Pan setting
         xo=length(Y)/2; % Initial Pan setting
         dx=length(Y)/4; % Initial Zoom setting
         AUTOZERO=0;
      else % if not isscalar
        % Must be separate x and y data vectors
        X=DataMatrix;
        Y=PeakD;
        PeakDensity=20;   
        % Estimate approx number of points in a peak half-width
        WidthPoints=length(Y)/PeakDensity;  
        SlopeThreshold=WidthPoints^-2;  
        AmpThreshold=abs(min(Y)+0.1*(max(Y)-min(Y))); 
        SmoothWidth=round(WidthPoints/3);  
        FitWidth=round(WidthPoints/3); 
        if FitWidth>100,FitWidth=100;end   % Keep FitWidth below 100
        if SmoothWidth>100,SmoothWidth=100;end   % Keep SmoothWidth below 100
        xo=length(Y)/2; % Initial Pan setting
        dx=length(Y)/4; % Initial Zoom setting
      end  % if isscalar
      AUTOZERO=0;
    case 3
      % Must be separate x and y data vectors plus a peak density estimate. 
      X=DataMatrix;
      Y=PeakD;
      % Calculate values of peak detection parameters
      % arguments based on the peak density, PeakD
      PeakDensity=AmpT;    
      % Estimate approx number of points in a peak half-width
      WidthPoints=length(Y)/PeakDensity;  
      SlopeThreshold=WidthPoints^-2;  
      AmpThreshold=abs(min(Y)+0.02*(max(Y)-min(Y)));  
      SmoothWidth=round(WidthPoints/3);  
      FitWidth=round(WidthPoints/3); 
      if FitWidth>100,FitWidth=100;end   % Keep FitWidth below 100
      if SmoothWidth>100,SmoothWidth=100;end   % Keep SmoothWidth below 100
      xo=length(Y)/2; % Initial Pan setting
      dx=length(Y)/4; % Initial Zoom setting
      AUTOZERO=0;
    case 6 
      % Must be one data matrix and all peak detection parameters 
      % specified in arguments 
      % If DataMatrix is in the wrong transposition, fix it.
      datasize=size(DataMatrix);
      if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
      X=DataMatrix(:,1); % Split matrix argument 
      Y=DataMatrix(:,2);
      SlopeThreshold=SlopeT;
      AmpThreshold=AmpT;
      SmoothWidth=SmoothW;
      FitWidth=FitW;
      xo=length(Y)/2; % Initial Pan setting
      dx=length(Y)/4; % Initial Zoom setting
      AUTOZERO=0;
    case 8
      % One data matrix, all peak detection parameters specified
      % in arguments, initial values of xcenter and xrange specified.
      % If DataMatrix is in the wrong transposition, fix it.
      datasize=size(DataMatrix);
      if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
      X=DataMatrix(:,1); % Split matrix argument 
      Y=DataMatrix(:,2);
      SlopeThreshold=SlopeT;
      AmpThreshold=AmpT;
      SmoothWidth=SmoothW;
      FitWidth=FitW;
      if xcenter<min(X),
          disp(['Lowest X value is ' num2str(min(X)) ]),
          xcenter=min(X)+xrange;
      end
      if xcenter>max(X),
          disp(['Highest X value is ' num2str(max(X)) ]),
          xcenter=max(X)-xrange;
      end
      xo=val2ind(X,xcenter);
      hirange=val2ind(X,xcenter+xrange./2);
      lorange=val2ind(X,xcenter-xrange./2);
      dx=(hirange-lorange);
      AUTOZERO=0;
    case 9
      % Like case 9, except initial AUTOZERO mode is specified
      % in arguments, initial values of xcenter and xrange specified.
      % If DataMatrix is in the wrong transposition, fix it.
      datasize=size(DataMatrix);
      if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
      X=DataMatrix(:,1); % Split matrix argument 
      Y=DataMatrix(:,2);
      SlopeThreshold=SlopeT;
      AmpThreshold=AmpT;
      SmoothWidth=SmoothW;
      FitWidth=FitW;
      xo=val2ind(X,xcenter);
      hirange=val2ind(X,xcenter+xrange./2);
      lorange=val2ind(X,xcenter-xrange./2);
      dx=(hirange-lorange);
      if xcenter<min(X),
          disp(['Lowest X value is ' num2str(min(X)) ]),
          xcenter=min(X)+xrange;
      end
      if xcenter>max(X),
          disp(['Highest X value is ' num2str(max(X)) ]),
          xcenter=max(X)-xrange;
      end
      AUTOZERO=MaxError;
    case 11 % last 3 options arguments provided
      datasize=size(DataMatrix);
      if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
      X=DataMatrix(:,1); % Split matrix argument 
      Y=DataMatrix(:,2);
      SlopeThreshold=SlopeT;
      AmpThreshold=AmpT;
      SmoothWidth=SmoothW;
      FitWidth=FitW;
      xo=val2ind(X,xcenter);
      hirange=val2ind(X,xcenter+xrange./2);
      lorange=val2ind(X,xcenter-xrange./2);
      dx=(hirange-lorange);
      if xcenter<min(X),
          disp(['Lowest X value is ' num2str(min(X)) ]),
          xcenter=min(X)+xrange;
      end
      if xcenter>max(X),
          disp(['Highest X value is ' num2str(max(X)) ]),
          xcenter=max(X)-xrange;
      end
        maxerror=MaxError;
        Positions=positions;
        Names=names;
        AUTOZERO=0;
    otherwise
      disp('Invalid number of arguments')
      disp('Expected forms are:')
      disp('ipeak(y);  % Data in single y vector')
      disp('ipeak(x,y);  % Data in separate x and y vectors')
      disp('ipeak(DataMatrix); % Data in two columns of DataMatrix')
      disp('ipeak(x,y,10), ipeak([x;y],10) or ipeak(y,10), specifying peak density')
      disp('ipeak(DataMatrix,0,.5,.0001,33,33);  specifying peak density, AmpT, SlopeT, SmoothW, FitW')
      disp('ipeak(DataMatrix,PeakD,AmpT,SlopeT,SmoothW,FitW,xcenter,xrange)')
      disp('ipeak(DataMatrix,PeakD,AmpT,SlopeT,SmoothW,FitW,xcenter,xrange,MaxError,positions,names)')
      beep
      return
end % switch nargin
% If necessary, flip the data vectors so that X increases
if X(1)>X(length(X)),
    disp('X-axis flipped.')
    X=fliplr(X);
    Y=fliplr(Y);
end
% ***********************************************
if FitWidth<2,FitWidth=2;end   % Keep FitWidth above 1 
PeakLabels=0; % Peak numbers only, no parameter labels, in upper window
PeakID=0; % Start with PeakID off
logplot=0; % Start with linear mode
plotcolor=0;  % Start with blue plot color
showpeak=1;  % Start with first peak under green cursor
valleymode=0;
% Plot the signal
P=findpeaks(X,Y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);
[xx,yy]=RedrawSignal(X,Y,xo,dx);
sizeP=size(P);
NumPeaks=sizeP(1);
P=MeasurePeaks(NumPeaks,X,Y,P,dx,SmoothWidth,FitWidth,AUTOZERO,valleymode);
% Attaches KeyPress test function to the figure.
set(gcf,'KeyPressFcn',@ReadKey)
uicontrol('Style','text')
% end of outer function
% ----------------------------SUBFUNCTIONS--------------------------------
function ReadKey(obj,eventdata)
% Interprets key presses from the Figure window. When a key is pressed, 
% executes the code in the corresponding section in the SWITCH statement.
% Note: If you don't like my key assignments, you can change the numbers
% in the case statements here to re-assign that function to any other key.
% If you press a key that has not yet been assigned to a function, it
% displays the key code number in the command window so you can easily
% add that to the SWITCH statement to add your own custom functions.
global X Y xx yy xo dx SlopeThreshold AmpThreshold SmoothWidth FitWidth plotcolor
global PeakLabels PeakID Names Positions maxerror SavedSignal oldAmpThreshold
global logplot P AUTOZERO showpeak valleymode
key=get(gcf,'CurrentCharacter');
if isscalar(key),
  ly=length(Y);
  switch double(key),
    case 29
        % Pans down when left arrow pressed.
        xo=xo+dx/10;
        if xo>ly,xo=ly;end
        [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 28
        % Pans up when right arrow pressed. 
         xo=xo-dx/10;
         if xo<1,xo=1;end
         [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 91
        % Nudge down 1 point when [ pressed.
         xo=xo-1;
         if xo<1,xo=1;end
         [xx,yy]=RedrawSignal(X,Y,xo,dx); 
    case 93
        % Nudge up 1 point when ] pressed.    
        xo=xo+1;
        if xo>ly,xo=ly;end
        [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 46
        % Pans down when < key pressed.
        xo=xo+dx/2;
        if xo>ly,xo=ly;end
        [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 44
        % Pans up when > key pressed.
        xo=xo-dx/2;
        if xo<1,xo=1;end
        [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 31
        % Zooms out when up arrow pressed.
        dx=dx+dx/10;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 30
        % Zooms in when down arrow pressed.
        dx=dx-dx/10;
        if dx<2,dx=2;end
        [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 47
        % Zooms out when / pressed.
        dx=dx+ly/50;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 39
        % Zooms in when ' pressed.
        dx=dx-ly/50;
        if dx<2,dx=2;end
        [xx,yy]=RedrawSignal(X,Y,xo,dx);      
     case 27 % When 'ESC' key is pressed, resets pan and zoom
         xo=length(Y)/2; % Initial Pan setting
         dx=length(Y)/4; % Initial Zoom setting
         [xx,yy]=RedrawSignal(X,Y,xo,dx);
     case 13  % Change plot color when Return (Enter) key pressed
         plotcolor=plotcolor+1;
         if plotcolor==6, plotcolor=0;end
         [xx,yy]=RedrawSignal(X,Y,xo,dx);
     case 98
        % When 'b' key is pressed, user clicks graph 
        % to enter background points, then graph re-drawn.
        SavedSignal=Y;
        oldAmpThreshold=AmpThreshold;
        BaselinePoints=input('Number of baseline points to click: ');
        if isempty(BaselinePoints),BaselinePoints=8;end
        AmpThreshold=input('Amplitude Threshold: ');
        if isempty(AmpThreshold),AmpThreshold=oldAmpThreshold;end
        % Acquire background points from user mouse clicks
        subplot(2,1,2)
        title(['Click on ' num2str(BaselinePoints) ' points on the baseline between the peaks.'])
        bX=[];bY=[];
        for g=1:BaselinePoints;
           [clickX,clickY] = ginput(1);
           bX(g)=clickX;
           bY(g)=clickY;
           xlabel(['Baseline point '  num2str(g) ' / ' num2str(BaselinePoints) ])
        end
        yy=Y;
        for k=1:length(bX)-1,
           fp=val2ind(X,bX(k)); % First point in segment
           lp=val2ind(X,bX(k+1));  % Last point in segment
           % Subtract piecewise linear background from Y
           yy(fp:lp)=Y(fp:lp)-((bY(k+1)-bY(k))/(bX(k+1)-bX(k))*(X(fp:lp)-bX(k))+bY(k));
        end
        Y=yy;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);  
    case 103
          % When 'g' key is pressed, restores signal background and AmpThreshold. 
          Y=SavedSignal;
          AmpThreshold=oldAmpThreshold;
          [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 97
        % When 'a' key is pressed, increases "AmpThreshold" by 10%
        AmpThreshold=AmpThreshold+.1*AmpThreshold;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);      
    case 122
        % When 'z' key is pressed, decreases "AmpThreshold" by 10%
        AmpThreshold=AmpThreshold-.1*AmpThreshold;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);      
    case 115 % When 's' key is pressed, increases "SlopeThreshold" by 10%
         SlopeThreshold=SlopeThreshold+.1*SlopeThreshold;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);   
    case 120 % When 'x' key is pressed, decreases "SlopeThreshold" by 10%
         SlopeThreshold=SlopeThreshold-.1*SlopeThreshold;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);   
    case 100
        % When 'd' key is pressed, increases "SmoothWidth" by 1 or 10%
        if SmoothWidth>20,
            SmoothWidth=round(SmoothWidth+.1.*SmoothWidth);
        else
            SmoothWidth=SmoothWidth+1;
        end
        [xx,yy]=RedrawSignal(X,Y,xo,dx);   
    case 99
        % When 'c' key is pressed, decreases "SmoothWidth" by 1 or 10%
        if SmoothWidth>20,
            SmoothWidth=round(SmoothWidth-.1.*SmoothWidth);
        else
            SmoothWidth=SmoothWidth-1;
        end
        if SmoothWidth<1, SmoothWidth=1;end
        [xx,yy]=RedrawSignal(X,Y,xo,dx);   
    case 102
        % When 'f' key is pressed, increases "FitWidth" by 1 or 10%
        if FitWidth>20,
            FitWidth=round(FitWidth+.1.*FitWidth);
        else
            FitWidth=FitWidth+1;
        end
        [xx,yy]=RedrawSignal(X,Y,xo,dx);   
    case 118
        % When 'v' key is pressed, decreases "FitWidth" by 1 or 10%
        if FitWidth>20,
            FitWidth=round(FitWidth-.1.*FitWidth);
        else
            FitWidth=FitWidth-1;
        end
        % **************************************************
         if FitWidth<2, FitWidth=2;end
        [xx,yy]=RedrawSignal(X,Y,xo,dx); 
    case 45
         % When '-' key is pressed, invert the signal
          Y=-Y;
          [xx,yy]=RedrawSignal(X,Y,xo,dx); 
          disp('Signal was inverted.')
    case 48
    % When '0' (zero) key is pressed, subtracts minimum from entire signal
    % (to remove positive or negative offset).
          Y=Y-min(Y);
          [xx,yy]=RedrawSignal(X,Y,xo,dx); 
          disp('Mininum signal set to zero.')
    case 114
        % When 'r' key is pressed, prints a report listing current 
        % settings and peak table.
        disp('--------------------------------------------------------')
        disp(['Amplitude Threshold (AmpT) = ' num2str(AmpThreshold) ] )
        disp(['Slope Threshold (SlopeT) = ' num2str(SlopeThreshold) ] )
        disp(['Smooth Width (SmoothW) = ' num2str(SmoothWidth) ' points' ] )
        disp(['Fit Width (FitW) = ' num2str(FitWidth) ' points' ] )
        sizeP=size(P);
        NumPeaks=sizeP(1);
        window=max(xx)-min(xx);
        if AUTOZERO,
            disp('Autozero ON')
            disp([ 'Window span: ' num2str(window) ]);
        else
            disp('Autozero OFF')
        end
          if valleymode,
              disp('        Valley#     Position     Height      Width          Area')
          else
              disp('          Peak#     Position     Height      Width          Area')
          end
        PP=MeasurePeaks(NumPeaks,X,Y,P,dx,SmoothWidth,FitWidth,AUTOZERO,valleymode);
        disp(PP)
     case 112
          % When 'p' key is pressed, prints out peak table
          disp('--------------------------------------------------------')
          sizeP=size(P);
          NumPeaks=sizeP(1);
          window=max(xx)-min(xx);
          if AUTOZERO,
              disp('Autozero ON')
              disp([ 'Window span: ' num2str(window) ' units'])
          else
              disp('Autozero OFF')
          end
          if valleymode,
              disp('        Valley#     Position     Height      Width          Area')
          else
              disp('          Peak#     Position     Height      Width          Area')
          end
          PP=MeasurePeaks(NumPeaks,X,Y,P,dx,SmoothWidth,FitWidth,AUTOZERO,valleymode);
          disp(PP)
    case 107
        % When 'k' key is pressed, prints out table of keyboard commands
        disp('KEYBOARD CONTROLS:')
        disp(' Pan left and right..........Coarse pan: < and >')   
        disp('                             Fine pan: left and right cursor arrows')
        disp('                             Nudge: [ ] ')
        disp(' Zoom in and out.............Coarse zoom: / and "  ') 
        disp('                             Fine zoom: up and down cursor arrows')
        disp(' Resets pan and zoom.........ESC')
        disp(' Change plot color...........Enter  (cycles through standard colors)')
        disp(' Adjust AmpThreshold.........A,Z (Larger values ignore short peaks)')
        disp(' Adjust SlopeThreshold.......S,X (Larger values ignore broad peaks)')
        disp(' Adjust SmoothWidth..........D,C (Larger values ignore sharp peaks}')
        disp(' Adjust FitWidth.............F,V (Adjust to cover just top part of peaks')
        disp(' Baseline correction:        B, enter # points, then click baseline ')
        disp(' Restore original signal.....G  to cancel previous background subtraction')
        disp(' Invert signal...............-  Invert (negate) the signal (flip + and -)')
        disp(' Set minimum to zero.........0  (Zero) Sets minumun signal to zero') 
        disp(' Toggle log y mode OFF/ON....Y  Plot log Y axis on lower graph')        
        disp(' Toggle autozero OFF/ON......T  Auto background subtraction on upper graph')  
        disp(' Toggle valley mode OFF/ON...U  Switch to valley mode')          
        disp(' Print report................R  prints Peak table and parameters')        
        disp(' Step through peaks..........Space/Tab  Jumps to next/previous peak')
        disp(' Print peak table............P  Peak #, Position, Height, Width, Area')
        disp(' Normal peak fit.............N  Fit peaks in upper window with peakfit.m')
        disp(' Multiple peak fit...........M  Fit all peaks in signal with peakfit.m')
        disp(' Print keyboard commands.....K  prints this list')
        disp(' Print findpeaks arguments...Q  prints findpeaks function with arguments')
        disp(' Print ipeak arguments.......W  prints ipeak function with all arguments')    
        disp(' Peak labels ON/OFF..........L  displays peak parameters on upper graph')
        disp(' Peak ID ON/OFF..............I  Identifies closest peaks in Names database.')
        disp(' Print table of peak IDs.....O  Prints Name, Position, Error, Amplitude')
    case 113
        % When 'Q' is pressed, prints findpeaks function with arguments    
        disp(['findpeaks(x,y,'  num2str(SlopeThreshold) ',' num2str(AmpThreshold)  ',' num2str(SmoothWidth)  ',' num2str(FitWidth) ',3)'] )
    case 119
        % When 'W' is pressed, prints ipeak function with arguments   
          center=(max(xx)+min(xx))/2;
          window=max(xx)-min(xx); 
        disp(['ipeak(DataMatrix,0,'  num2str(AmpThreshold) ',' num2str(SlopeThreshold)  ',' num2str(SmoothWidth)  ',' num2str(FitWidth) ',' num2str(center) ',' num2str(window) ')'] )
    case 105
        % When 'I' is pressed, toggles on/off PeakID in upper panel
        if PeakID==0,
            PeakID=1;
            % load DataTable
            % disp([ 'Loaded "DataTable" from disk. Number of Names:' num2str(length(Positions)) ] )
            % disp(['Position range: ' num2str(min(Positions)) '-' num2str(max(Positions)) ] )
        else
            PeakID=0;
        end
        [xx,yy]=RedrawSignal(X,Y,xo,dx); 
    case 108
         % When 'L' is pressed, toggles on/off peak labels in upper panel
        if PeakLabels==0,
            PeakLabels=1;
        else
            PeakLabels=0;
        end
        [xx,yy]=RedrawSignal(X,Y,xo,dx); 
    case 121
         % When 'Y' is pressed, toggles on/off log plot mode
        if logplot==0,
            logplot=1;
        else
            logplot=0;
        end
        [xx,yy]=RedrawSignal(X,Y,xo,dx); 
     case 117
         % When 'U' is pressed, toggles valleymode on/off 
        if valleymode==0,
            valleymode=1;
        else
            valleymode=0;
        end
        [xx,yy]=RedrawSignal(X,Y,xo,dx); 
    case 111
          % When 'o' is pressed, prints table of identified peaks
          if PeakID,
              disp('      Name          Position      Error         Amplitude')  %  Print out column lables for table
              for n=1:length(P(:,2)),
                  % m=index of the cloest match in Positions
                  m=val2ind(Positions,P(n,2));
                  % Error=difference between detected peak and nearest
                  % peak in table
                  Error=abs(P(n,2)-Positions(m));
                  if Error<maxerror, % Only identify the peaks if the error is less than MaxError
                      disp([Names(m) Positions(m) Error P(n,3)]); % Print out one line of Positions and Errors table
                  end % if error
              end  % for n
          end  % if PeakID

    case 110
         % When 'N' is pressed, applies peakfit function only to peaks in
         % the upper window (up to 6 peaks).
         % [xx,yy]=RedrawSignal(X,Y,xo,dx); 
         sizeP=size(P);
         NumPeaksUW=sizeP(1);
         if NumPeaksUW>1,
             PUW=[]; % PUW=table of peaks in upper window
             for peak=1:NumPeaksUW, % NumPeaksUW=number of peaks in upper window
                 if P(peak,2)>min(xx),
                     if P(peak,2)<max(xx),
                         PUW=[PUW;P(peak,:)];
                     end %  if P(peak,2)<max(xx),
                 end %  if P(peak,2)>min(xx),
             end % for peak=1:length(P),
         else
             PUW=P;
         end
          sizePUW=size(PUW);
          NumPeaksUW=sizePUW(1);
          center=(max(xx)+min(xx))/2;
          window=max(xx)-min(xx); 
          extra=1;
          disp('1=Gaussian (default), 2=Lorentzian, 3=logistic, 4=Pearson'); 
          disp('5=exponentionally broadened Gaussian, 6=equal-width Gaussians');
          disp('7=Equal-width Lorentzians, 8=exponentionally broadened equal-width Gaussian');
          Shape=input('Peak shape (1-8): ');
          if isempty(Shape),Shape=1;end
          NumTrials=input('Number of trials: ');
          if isempty(NumTrials),NumTrials=1;end
          if Shape==4||Shape==5||Shape==8,
             extra=input('Extra parameter: ');
          end % if Shape==4||Shape==5||Shape==8,
          if NumTrials>1,disp(['Best of ' num2str(NumTrials) ' trial fits.' ]), end       
          switch NumPeaksUW
            case 1
              startvector=[PUW(1,2) PUW(1,4)];
            case 2
              startvector=[PUW(1,2) PUW(1,4) PUW(2,2) PUW(2,4)];    
            case 3
              startvector=[PUW(1,2) PUW(1,4) PUW(2,2) PUW(2,4) PUW(3,2) PUW(3,4)];    
            case 4
              startvector=[PUW(1,2) PUW(1,4) PUW(2,2) PUW(2,4) PUW(3,2) PUW(3,4) PUW(4,2) PUW(4,4)];
            case 5
              startvector=[PUW(1,2) PUW(1,4) PUW(2,2) PUW(2,4) PUW(3,2) PUW(3,4) PUW(4,2) PUW(4,4) PUW(5,2) PUW(5,4)];        
            case 6
              startvector=[PUW(1,2) PUW(1,4) PUW(2,2) PUW(2,4) PUW(3,2) PUW(3,4) PUW(4,2) PUW(4,4) PUW(5,2) PUW(5,4) PUW(6,2) PUW(6,4)];        
          end % switch NumPeaksUW
          [FitResults,MeanFitError]=peakfit([xx,yy],center,window,NumPeaksUW,Shape,extra,NumTrials,startvector,AUTOZERO);
          switch Shape
              case 1
                  ShapeString='Gaussian';
              case 2
                  ShapeString='Lorentzian';
              case 3
                  ShapeString='logistic';
              case 4
                  ShapeString='Pearson7';
              case 5
                  ShapeString='ExpGaussian';
              case 6
                  ShapeString='Equal width Gaussians';
              case 7
                  ShapeString='Equal width Lorentzians';
              case 8
                  ShapeString='Equal-width ExpGauss.';
              otherwise
                  ShapeString='';
          end % switch Shape
        disp(['Least-squares fit to ' ShapeString ' peak model' ])
          disp(['Fitting Error ' num2str(MeanFitError) '%'])
        disp('          Peak#     Position     Height      Width         Area  ') 
        for peak=1:NumPeaksUW,FitResults(peak,1)=PUW(peak,1);end
        disp(FitResults(:,1:5))
     case 116
        % When 't' key is pressed, toggles AUTOZERO mode
        if AUTOZERO,
           AUTOZERO=0;
           [xx,yy]=RedrawSignal(X,Y,xo,dx); 
        else
           AUTOZERO=1;
           [xx,yy]=RedrawSignal(X,Y,xo,dx); 
        end
      case 101
          % When 'e' is pressed,
          disp('-------------------------------------------------------------')
    case 109
        % When 'M' is pressed, applies peakfit function to all peaks
        AllFitResults=[];
        [xx,yy]=RedrawSignal(X,Y,xo,dx);
        sizeP=size(P);
        NumPeaks=sizeP(1);
        if NumPeaks>2,
            center=(max(xx)+min(xx))/2;
            window=max(xx)-min(xx);
            extra=1;
            disp('1=Gaussian (default), 2=Lorentzian, 3=logistic, 4=Pearson');
            disp('5=exponentionally broadened Gaussian, 6=equal-width Gaussians');
            disp('7=Equal-width Lorentzians, 8=exponentionally broadened equal-width Gaussian');
            Shape=input('Peak shape (1-8): ');
            if isempty(Shape),Shape=1;end
            NumTrials=input('Number of trials (1-100): ');
            if isempty(NumTrials),NumTrials=1;end
            if Shape==4||Shape==5||Shape==8,
                extra=input('Extra parameter: ');
            end
            switch Shape
                case 1
                    ShapeString='Gaussian';
                case 2
                    ShapeString='Lorentzian';
                case 3
                    ShapeString='logistic';
                case 4
                    ShapeString='Pearson7';
                case 5
                    ShapeString='ExpGaussian';
                case 6
                    ShapeString='Equal width Gaussians';
                case 7
                    ShapeString='Equal width Lorentzians';
                case 8
                    ShapeString='Equal-width ExpGauss.';
                otherwise
                    ShapeString='';
            end % switch Shape
            disp(['Multiple Least-squares fit to ' ShapeString ' peak model' ]);
            disp('          Peak#     Position     Height      Width         Area          Error')
            for peak=1:NumPeaks-1,
                xcenter=P(peak,2);
                xrange=8*P(peak,4);
                xo=val2ind(X,xcenter);
                hirange=val2ind(X,xcenter+xrange./2);
                lorange=val2ind(X,xcenter-xrange./2);
                dx=(hirange-lorange);
                [xx,yy]=RedrawSignal(X,Y,xo,dx);
                PP=findpeaks(xx,yy,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);
                sizeP=size(PP);
                NumPeaksUW=sizeP(1);% Number of peaks in Upper Window
                switch NumPeaksUW
                    case 1
                        startvector=[PP(1,2) PP(1,4)];
                    case 2
                        startvector=[PP(1,2) PP(1,4) PP(2,2) PP(2,4)];
                    case 3
                        startvector=[PP(1,2) PP(1,4) PP(2,2) PP(2,4) PP(3,2) PP(3,4)];
                    case 4
                        startvector=[PP(1,2) PP(1,4) PP(2,2) PP(2,4) PP(3,2) PP(3,4) PP(4,2) PP(4,4)];
                    case 5
                        startvector=[PP(1,2) PP(1,4) PP(2,2) PP(2,4) PP(3,2) PP(3,4) PP(4,2) PP(4,4) PP(5,2) PP(5,4)];
                    case 6
                        startvector=[PP(1,2) PP(1,4) PP(2,2) PP(2,4) PP(3,2) PP(3,4) PP(4,2) PP(4,4) PP(5,2) PP(5,4) PP(6,2) PP(6,4)];
                end % switch NumPeaksUW
                
                [FitResults,MeanFitError]=peakfit([xx,yy],center,window,NumPeaksUW,Shape,extra,NumTrials,startvector,0);
                % for peak=1:NumPeaks,FitResults(peak,1)=PUW(peak,1);end
                for fittrial=1:NumPeaksUW,  % Number of peaks in Upper Window
                    AllFitResults=[AllFitResults;[val2ind(P(:,2),FitResults(fittrial,2)) FitResults(fittrial,2:5) MeanFitError]];
                end % for fittrial
            end % peak=1:NumPeaks-1,
            % Select the best for for each peak
            SortedResults=sortrows(AllFitResults,2); % Sort AllFitResults by position (column 2)
            SizeAllFitResults=size(AllFitResults);
            NumFits=SizeAllFitResults(1);
            BestFits=[];
            for FitNumber=1:NumPeaks,
                FirstColumn=min(val2ind(SortedResults(:,1),FitNumber));
                LastColumn=max(val2ind(SortedResults(:,1),FitNumber));
                SelectedSection=SortedResults(FirstColumn:LastColumn,:);
                SortedSection=sortrows(SelectedSection,6);
                BestRow=SortedSection(1,:);
                BestFits=[BestFits;BestRow];
            end
            disp(BestFits)
        else disp('Too few peaks detected; use the Normal curve fit instead.')
        end % if  NumPeaks>2,
      case 32
          % When Spacebar is pressed, jumps to next peak
          if valleymode,
              P=findvalleys(X,Y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth);
          else
              P=findpeaks(X,Y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);
          end
          sizeP=size(P);
          NumPeaks=sizeP(1);
          showpeak=showpeak+1;
          if showpeak>NumPeaks,showpeak=1;end
          center=P(showpeak,2);
          xo=val2ind(X,center);
          [xx,yy]=RedrawSignal(X,Y,xo,dx);
      case 9
          % When Tab is pressed, jumps to previous peak
          if valleymode,
              P=findvalleys(X,Y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);
          else
              P=findpeaks(X,Y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);
          end
          sizeP=size(P);
          NumPeaks=sizeP(1);
          showpeak=showpeak-1;
          if showpeak>NumPeaks,showpeak=1;end
          if showpeak<1,showpeak=NumPeaks;end
          center=P(showpeak,2);
          xo=val2ind(X,center);
          [xx,yy]=RedrawSignal(X,Y,xo,dx);
      otherwise  
         UnassignedKey=double(key)
       disp('Press k to print out list of keyboard commands')
   end % switch double(key),
end % if ischar(key),
% ----------------------------------------------------------------------    
function [xx,yy]=RedrawSignal(X,Y,xo,dx)
% Plots the entire signal (X,Y) in the lower half of the plot window and an
% isolated segment (xx,yy) in the upper half, controlled by Pan and Zoom
% keys.
global SlopeThreshold AmpThreshold SmoothWidth FitWidth PeakLabels valleymode
global PeakID Names Positions maxerror P plotcolor logplot AUTOZERO 
Startx=round(xo-(dx/2));
Endx=abs(round(xo+(dx/2)-1));
if (Endx-Startx)<SmoothWidth,Endx=Startx+SmoothWidth;end
if Endx>length(Y),Endx=length(Y);end
if Startx<1,Startx=1;end
PlotRange=Startx:Endx;
if (Endx-Startx)<5, PlotRange=xo:xo+5;end
xx=X(PlotRange);
yy=Y(PlotRange);
hold off
% clf
% Plots isolated segment (xx,yy) in the upper half
switch plotcolor
    case 0
        color='b.';
    case 1
        color='g.';
    case 2
        color='r.';
    case 3
        color='c.';
    case 4
        color='m.';
    case 5
        color='k.';
end
% auto-zero operation
if AUTOZERO==1,
    X1=min(xx);
    X2=max(xx);
    Y1=mean(yy(1:length(xx)/20));
    Y2=mean(yy((length(xx)-length(xx)/20):length(xx)));
    yy=yy-((Y2-Y1)/(X2-X1)*(xx-X1)+Y1);
end % if

figure(1);
hold off
if logplot,
    semilogy(xx,abs(yy),color)  % Graph the signal with linear Y axis
else
    subplot(2,1,1);plot(xx,yy,color);  % Graph the signal with linear Y axis
end
% subplot(2,1,1);plot(xx,yy,color);

hold off
if AUTOZERO==1
    if valleymode
         title('iPeak 3.9.   Valley mode. Autozero ON.  Press K for keyboard commands')
    else
         title('iPeak 3.9.   Peak mode.   Autozero ON.  Press K for keyboard commands')
    end
else
    if valleymode
         title('iPeak 3.9.   Valley mode. Autozero OFF.  Press K for keyboard commands')
    else
         title('iPeak 3.9.   Peak mode.   Autozero OFF.  Press K for keyboard commands')
    end
end

axis([X(Startx(1)) X(Endx(1)) min(yy) max(yy)+(max(yy)-min(yy))/10]);
xlabel('Space/Tab: next/previous peak.  Mode: U  Autozero: T   Log/linear: Y  Report: R')

% Bottom half of the figure shows full signal
subplot(2,1,2);cla
switch plotcolor
    case 0; color='b';
    case 1; color='g';
    case 2; color='r';
    case 3; color='c';
    case 4; color='m';
    case 5; color='k';
end
hold off
if logplot,
    semilogy(X,abs(Y),color)  % Graph the signal with linear Y axis
else
    plot(X,Y,color)  % Graph the signal with linear Y axis
end
if valleymode,
  P=findvalleys(X,Y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);
else
  P=findpeaks(X,Y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);   
end
title('AmpT: A/Z   SlopeT: S/X   SmoothW: D/C    FitW: F/V   Background: B'  )
if logplot,
    ylabel('Log y mode')
    xlabel(['    AmpT: ' num2str(AmpThreshold) '     SlopeT: ' num2str(SlopeThreshold) '    SmoothW: ' num2str(SmoothWidth) '    FitW: ' num2str(FitWidth) ])
    axis([X(1) X(length(X)) min(abs(Y)) max(Y)]); % Update plot
else
    ylabel('Linear mode')
    xlabel(['    AmpT: ' num2str(AmpThreshold) '     SlopeT: ' num2str(SlopeThreshold) '    SmoothW: ' num2str(SmoothWidth) '    FitW: ' num2str(FitWidth) ])
    axis([X(1) X(length(X)) min(Y) max(Y)]); % Update plot
end
text(P(:,2),P(:,3),num2str(P(:,1)))  % Number the peaks found on the lower graph
hold on
% Mark the zoom range on the full signal with two magenta dotted vertical lines
center=X(round(xo));
checkzero=abs(Y);
checkzero(~checkzero)=NaN; % Find smallest non-zero value
plot([min(xx) min(xx)],[min(checkzero) max(Y)],'m--')
plot([max(xx) max(xx)],[min(checkzero) max(Y)],'m--')
plot([center center],[min(checkzero) max(Y)],'g-')
% Number the peaks found on the upper graph
subplot(2,1,1);
if valleymode,
  PP=findvalleys(xx,yy,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);
else
  PP=findpeaks(xx,yy,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);   
end
if PeakLabels,
    % Label the peaks on the upper graph with number, position, height, and
    % width
    % auto-zero operation
    if AUTOZERO==1,
        X1=min(xx);
        X2=max(xx);
        Y1=mean(yy(1:length(xx)/20));
        Y2=mean(yy((length(xx)-length(xx)/20):length(xx)));
        yy=yy-((Y2-Y1)/(X2-X1)*(xx-X1)+Y1);
    end % if AUTOZERO
    topaxis=axis;
    yrange=topaxis(4)-topaxis(3);
    pos1=.1*yrange;
    pos2=.2*yrange;
    pos3=.3*yrange;
    pos4=.4*yrange;
    text(P(:,2),P(:,3)-pos1,num2str(P(:,1)))
    text(PP(:,2),PP(:,3)-pos2,num2str(PP(:,2)))
    text(PP(:,2),PP(:,3)-pos3,num2str(PP(:,3)))
    text(PP(:,2),PP(:,3)-pos4,num2str(PP(:,4)))
else
    topaxis=axis;
    yrange=topaxis(4)-topaxis(3);
    pos1=.1*yrange;
    % Number the peaks on the upper graph
    text(P(:,2),P(:,3)-pos1,num2str(P(:,1)))
end
hold on
lyy=min(yy);
uyy=max(yy)+(max(yy)-min(yy))/10;
if lyy<uyy;
    axis([X(Startx(1)) X(Endx(1)) lyy uyy ]);
end
center=X(round(xo));
hold on;plot([center center],[lyy uyy],'g-')
% Draw red best-fit line through peak tops in upper windows.
if PP(1)>0, % if any peaks are detected
    sizePP=size(PP);
    lengthPP=sizePP(1);
    for PeakNumber=1:lengthPP,
        subplot(2,1,1);
        if PeakNumber>lengthPP,PeakNumber=lengthPP;end
        n1=round(val2ind(xx,PP(PeakNumber,2))-FitWidth/2);
        n2=round(val2ind(xx,PP(PeakNumber,2))+FitWidth/2);
        if n1<1, n1=1;end
        if n2>length(yy), n2=length(yy);end
        PlotRange=n1:n2;
        xxx=xx(PlotRange);
        yyy=yy(PlotRange);
        if valleymode,
            [coef,S]=polyfit(xxx,yyy,2);  % Fit parabola to sub-group with centering and scaling
            c1=coef(3);c2=coef(2);c3=coef(1);
            subplot(2,1,1);
            plotspace=linspace(min(xxx),max(xxx));
            plot(plotspace,c3*plotspace.^2+c2*plotspace+c1,'r');
        else
            % Fit parabola to log10 of sub-group
            [coef,S,MU]=polyfit(xxx,log(abs(yyy)),2);
            c1=coef(3);c2=coef(2);c3=coef(1);
            % Compute peak position and height of fitted parabola
            PeakX=-((MU(2).*c2/(2*c3))-MU(1));
            PeakY=exp(c1-c3*(c2/(2*c3))^2);
            MeasuredWidth=norm(MU(2).*2.35482/(sqrt(2)*sqrt(-1*c3)));
            subplot(2,1,1);
            plotspace=linspace(min(xxx),max(xxx));
            plot(plotspace,PeakY.*gaussian(plotspace,PeakX,MeasuredWidth),'r');
        end
    end   % for PeakNumber
    % Place a label in the upper left corner with peak number, position,
    % height, and width of the peak closest to the center of the window.
    PeakAtCenter=val2ind(P(:,2),center);
    % auto-zero operation
    if AUTOZERO==1,
        X1=min(xx);
        X2=max(xx);
        Y1=mean(yy(1:length(xx)/20));
        Y2=mean(yy((length(xx)-length(xx)/20):length(xx)));
        yy=yy-((Y2-Y1)/(X2-X1)*(xx-X1)+Y1);
    end % if AUTOZERO
    n1=round(val2ind(xx,P(PeakAtCenter,2))-FitWidth/2);
    n2=round(val2ind(xx,P(PeakAtCenter,2))+FitWidth/2);
    if n1<1, n1=1;end
    if n2>length(yy), n2=length(yy);end
    FitRange=n1:n2;
    xxx=xx(FitRange);
    yyy=yy(FitRange);
    if valleymode,
        [coef]=polyfit(xxx,yyy,2);  % Fit parabola to sub-group with centering and scaling
        c1=coef(3);c2=coef(2);c3=coef(1);
        PeakX=-c2/(2*c3);
        PeakY=(c1-(c2*c2/(4*c3)));
        MeasuredWidth=norm(2.35482/(sqrt(2)*sqrt(-1*c3)));
    else
        % Fit parabola to log10 of sub-group
        [coef,S,MU]=polyfit(xxx,log(abs(yyy)),2);
        c1=coef(3);c2=coef(2);c3=coef(1);
        % Compute peak position and height or fitted parabola
        PeakX=-((MU(2).*c2/(2*c3))-MU(1));
        PeakY=exp(c1-c3*(c2/(2*c3))^2);
        MeasuredWidth=norm(MU(2).*2.35482/(sqrt(2)*sqrt(-1*c3)));
    end
    startx=min(xx)+(max(xx)-min(xx))./20;
    dyy=(max(yy)-min(yy))./10;
    starty=max(yy)-dyy;
    if valleymode,
        text(startx,starty+dyy/2,['Valley ' num2str(PeakAtCenter) ] );
    else
        text(startx,starty+dyy/2,['Peak ' num2str(PeakAtCenter) ] );
    end
    topaxis=axis;
    yrange=topaxis(4)-topaxis(3);
    pos1=.1*yrange;
    pos2=.2*yrange;
    pos3=.3*yrange;
    text(startx,starty+dyy/2-pos1,['Position: ' num2str(PeakX)])
    text(startx,starty+dyy/2-pos2,['Height: ' num2str(PeakY)])
    text(startx,starty+dyy/2-pos3,['Width: ' num2str(MeasuredWidth)])
    % Add peak identification if peak identification mode is ON and
    % information provided in arguments 9, 10, and 11.
    if PeakID, % If peak identification mode is ON
        for n=1:length(PP(:,2)),
            %      [PP(n,2) Positionsv(n)]
            m=val2ind(Positions,PP(n,2)); % m=index of the cloest match in Positions
            Error=abs(PP(n,2)-Positions(m)); % Error=difference between detected peak and nearest peak in table
            if Error<maxerror, % Only identify the peaks if the error is less than MaxError
                text(PP(n,2),PP(n,3),Names(m)); % Label the graph peaks with element names
            end % if error
        end  % for n
    end  % if PeakID
end  % if any peaks are detected
hold off
sizeP=size(P);
NumPeaks=sizeP(1);
P=MeasurePeaks(NumPeaks,X,Y,P,dx,SmoothWidth,FitWidth,AUTOZERO,valleymode);
%-----------------------------------------------------------------
function PP=MeasurePeaks(NumPeaks,X,Y,P,dx,SmoothWidth,FitWidth,AUTOZERO,valleymode)
PP=zeros(size(P));
for PeakNumber=1:NumPeaks,
    center=P(PeakNumber,2);
    xo=val2ind(X,center);
    Startx=round(xo-(dx/2));
    Endx=abs(round(xo+(dx/2)-1));
    if (Endx-Startx)<SmoothWidth,Endx=Startx+SmoothWidth;end
    if Endx>length(Y),Endx=length(Y);end
    if Startx<1,Startx=1;end
    PlotRange=Startx:Endx;
    if (Endx-Startx)<5, PlotRange=xo:xo+5;end
    xx=X(PlotRange);
    yy=Y(PlotRange);
    if AUTOZERO==1,
        X1=min(xx);
        X2=max(xx);
        Y1=mean(yy(1:length(xx)/20));
        Y2=mean(yy((length(xx)-length(xx)/20):length(xx)));
        yy=yy-((Y2-Y1)/(X2-X1)*(xx-X1)+Y1);
    end % if AUTOZERO
    n1=round(val2ind(xx,P(PeakNumber,2))-FitWidth/2);
    n2=round(val2ind(xx,P(PeakNumber,2))+FitWidth/2);
    if n1<1, n1=1;end
    if n2>length(yy), n2=length(yy);end
    FitRange=n1:n2;
    xxx=xx(FitRange);
    yyy=yy(FitRange);
    if valleymode,
        [coef]=polyfit(xxx,yyy,2);  % Fit parabola to sub-group xxx, yyy
        c1=coef(3);c2=coef(2);c3=coef(1);
        PeakX=-c2/(2*c3);
        PeakY=(c1-(c2*c2/(4*c3)));
        MeasuredWidth=norm(2.35482/(sqrt(2)*sqrt(-1*c3)));
        EstimatedArea=0;
    else
        [coef,S,MU]=polyfit(xxx,log(abs(yyy)),2);   % Fit parabola to log10 of sub-group
        c1=coef(3);c2=coef(2);c3=coef(1);
        % Compute peak position and height or fitted parabola
        PeakX=-((MU(2).*c2/(2*c3))-MU(1));
        PeakY=exp(c1-c3*(c2/(2*c3))^2);
        MeasuredWidth=norm(MU(2).*2.35482/(sqrt(2)*sqrt(-1*c3)));
        EstimatedArea=1.0646.*PeakY*MeasuredWidth;
    end
    PP(PeakNumber,:)=[PeakNumber PeakX PeakY MeasuredWidth EstimatedArea];
end   % for PeakNumber
% ----------------------------------------------------------------------
function P=findpeaks(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% function P=findpeaks(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% Function to locate the positive peaks in a noisy x-y time series data
% set.  Detects peaks by looking for downward zero-crossings
% in the first derivative that exceed SlopeThreshold.
% Returns list (P) containing peak number and position, 
% height, width, and area of each peak. Arguments "slopeThreshold",
% "ampThreshold" and "smoothwidth" control peak sensitivity.
% Higher values will neglect smaller features. "Smoothwidth" is
% the width of the smooth applied before peak detection; larger
% values ignore narrow peaks. "Peakgroup" is the number points 
% around the top part of the peak that are taken for measurement.
% The argument "smoothtype" determines the smooth algorithm:
%   If smoothtype=1, rectangular (sliding-average or boxcar) 
%   If smoothtype=2, triangular (2 passes of sliding-average)
%   If smoothtype=3, pseudo-Gaussian (3 passes of sliding-average)
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and 
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% T. C. O'Haver, 1995.  Version 4.1, Last revised September, 2011
if nargin~=7;smoothtype=1;end  % smoothtype=1 if not specified in argument
if smoothtype>3;smoothtype=3;end
if smoothtype<1;smoothtype=1;end 
smoothwidth=round(smoothwidth);
peakgroup=round(peakgroup);
d=fastsmooth(deriv(y),smoothwidth,smoothtype);
n=round(peakgroup/2+1);
P=[0 0 0 0 0];
vectorlength=length(y);
peak=1;
AmpTest=AmpThreshold;
for j=smoothwidth:length(y)-smoothwidth,
    if sign(d(j)) > sign (d(j+1)), % Detects zero-crossing
        if d(j)-d(j+1) > SlopeThreshold*y(j), % if slope of derivative is larger than SlopeThreshold
            if y(j) > AmpTest,  % if height of peak is larger than AmpThreshold
                xx=zeros(size(peakgroup));yy=zeros(size(peakgroup));
                for k=1:peakgroup, % Create sub-group of points near peak
                    groupindex=j+k-n+2;
                    if groupindex<1, groupindex=1;end
                    if groupindex>vectorlength, groupindex=vectorlength;end
                    xx(k)=x(groupindex);yy(k)=y(groupindex);
                end
                    [Height, Position, Width]=gaussfit(xx,yy);
                    PeakX=real(Position);   % Compute peak position and height of fitted parabola
                    PeakY=real(Height);
                    MeasuredWidth=real(Width);
                % if the peak is too narrow for least-squares technique to work
                % well, just use the max value of y in the sub-group of points near peak.
                if peakgroup<3,
                    PeakY=max(yy);
                    pindex=val2ind(yy,PeakY);
                    PeakX=xx(pindex(1));
                end
                % Construct matrix P. One row for each peak
                % detected, containing the peak number, peak
                % position (x-value) and peak height (y-value).
                P(peak,:) = [round(peak) PeakX PeakY MeasuredWidth  1.0646.*PeakY*MeasuredWidth];
                peak=peak+1;
  
            end
        end
    end
end
% ----------------------------------------------------------------------
function V=findvalleys(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% function P=findvalleys(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% Function to locate the valleys (mimnima) in a noisy x-y time series data
% set.  Detects valleys by looking for upward zero-crossings
% in the first derivative that exceed SlopeThreshold.
% Returns list (V) containing valley number and position, 
% depth, and width of each valley. Arguments "slopeThreshold",
% "ampThreshold" and "smoothwidth" control sensitivity.
% Higher values will neglect smaller features. "Smoothwidth" is
% the width of the smooth applied before valley detection; larger
% values ignore narrow features. "Peakgroup" is the number points 
% around the bottom part of the valley that are fit to a parabola to
% determine the valley vertex (x and y at lowest point) and width.
% The argument "smoothtype" determines the smooth algorithm:
%   If smoothtype=1, rectangular (sliding-average or boxcar) 
%   If smoothtype=2, triangular (2 passes of sliding-average)
%   If smoothtype=3, pseudo-Gaussian (3 passes of sliding-average)
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and 
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% T. C. O'Haver, Version 2.1, September, 2011
if nargin~=7;smoothtype=1;end  % smoothtype=1 if not specified in argument
if smoothtype>3;smoothtype=3;end
if smoothtype<1;smoothtype=1;end 
smoothwidth=round(smoothwidth);
peakgroup=round(peakgroup);
d=fastsmooth(deriv(y),smoothwidth,smoothtype);
n=round(peakgroup/2+1);
V=[0 0 0 0 0];
vectorlength=length(y);
peak=1;
AmpTest=AmpThreshold;
for j=smoothwidth:length(y)-smoothwidth,
    if sign(d(j)) < sign (d(j+1)), % Detects zero-crossing
        if d(j+1)-d(j) > SlopeThreshold*y(j), % if slope of derivative is larger than SlopeThreshold
            if y(j) > AmpTest,  % if height of valley is larger than AmpThreshold
                xx=zeros(size(peakgroup));yy=zeros(size(peakgroup));
                for k=1:peakgroup, % Create sub-group of points near valley
                    groupindex=j+k-n+1;
                    if groupindex<1, groupindex=1;end
                    if groupindex>vectorlength, groupindex=vectorlength;end
                    xx(k)=x(groupindex);yy(k)=y(groupindex);
                end
                [coef,S]=polyfit(xx,yy,2);  % Fit parabola to sub-group with centering and scaling
                c1=coef(3);c2=coef(2);c3=coef(1);
                valleyX=-c2/(2*c3);   % Compute valley position and height of fitted parabola
                valleyY=(c1-(c2*c2/(4*c3)));
                MeasuredWidth=norm(2.35482/(sqrt(2)*sqrt(-1*c3)));
                % if the valley is too narrow for least-squares technique to work
                % well, just use the min value of y in the sub-group of points near valley.
                if peakgroup<5,
                    valleyY=min(yy);
                    pindex=val2ind(yy,valleyY);
                    valleyX=xx(pindex(1));
                end
                % Construct matrix P. One row for each valley
                % detected, containing the valley number, valley
                % position (x-value) and valley depth (y-value).
                % Area is not measured for valleys, so put a zero
                V(peak,:) = [round(peak) valleyX valleyY MeasuredWidth 0];
                peak=peak+1;
            end
        end
    end
end
% ----------------------------------------------------------------------
function d=deriv(a)
% First derivative of vector using 2-point central difference.
%  T. C. O'Haver, 1988.
n=length(a);
d=zeros(size(a));
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1;
  d(j)=(a(j+1)-a(j-1)) ./ 2;
end
% ----------------------------------------------------------------------
function SmoothY=fastsmooth(Y,w,type,ends)
% fastsmooth(Y,w,type,ends) smooths vector Y with smooth 
%  of width w. Version 2.0, May 2008.
% The argument "type" determines the smooth type:
%   If type=1, rectangular (sliding-average or boxcar) 
%   If type=2, triangular (2 passes of sliding-average)
%   If type=3, pseudo-Gaussian (3 passes of sliding-average)
% The argument "ends" controls how the "ends" of the signal 
% (the first w/2 points and the last w/2 points) are handled.
%   If ends=0, the ends are zero.  (In this mode the elapsed 
%     time is independent of the smooth width). The fastest.
%   If ends=1, the ends are smoothed with progressively 
%     smaller smooths the closer to the end. (In this mode the  
%     elapsed time increases with increasing smooth widths).
% fastsmooth(Y,w,type) smooths with ends=0.
% fastsmooth(Y,w) smooths with type=1 and ends=0.
% Example:
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3)= [0 1 4 7 10 7 4 1 1 0]
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3,1,1)= [1 1 4 7 10 7 4 1 1 1]
%  T. C. O'Haver, May, 2008.
if nargin==2, ends=0; type=1; end
if nargin==3, ends=0; end
  switch type
    case 1
       SmoothY=sa(Y,w,ends);
    case 2   
       SmoothY=sa(sa(Y,w,ends),w,ends);
    case 3
       SmoothY=sa(sa(sa(Y,w,ends),w,ends),w,ends);
  end

function SmoothY=sa(Y,smoothwidth,ends)
w=round(smoothwidth);
SumPoints=sum(Y(1:w));
s=zeros(size(Y));
halfw=round(w/2);
L=length(Y);
for k=1:L-w,
   s(k+halfw-1)=SumPoints;
   SumPoints=SumPoints-Y(k);
   SumPoints=SumPoints+Y(k+w);
end
s(k+halfw)=sum(Y(L-w+1:L));
SmoothY=s./w;
% Taper the ends of the signal if ends=1.
  if ends==1,
    startpoint=(smoothwidth + 1)/2;
    SmoothY(1)=(Y(1)+Y(2))./2;
    for k=2:startpoint,
       SmoothY(k)=mean(Y(1:(2*k-1)));
       SmoothY(L-k+1)=mean(Y(L-2*k+2:L));
    end
    SmoothY(L)=(Y(L)+Y(L-1))./2;
  end
% ----------------------------------------------------------------------
  function [Height, Position, Width]=gaussfit(x,y)
% Converts y-axis to a log scale, fits a parabola
% (quadratic) to the (x,ln(y)) data, then calculates
% the position, width, and height of the
% Gaussian from the three coefficients of the
% quadratic fit.  This is accurate only if the data have
% no baseline offset (that is, trends to zero far off the
% peak) and if there are no zeros or negative values in y.
%
% Example 1: Simplest Gaussian data set
% [Height, Position, Width]=gaussfit([1 2 3],[1 2 1]) 
%    returns Height = 2, Position = 2, Width = 2
%
% Example 2: best fit to synthetic noisy Gaussian
% x=50:150;y=100.*gaussian(x,100,100)+10.*randn(size(x));
% [Height,Position,Width]=gaussfit(x,y) 
%   returns [Height,Position,Width] clustered around 100,100,100.
%
% Example 3: plots data set as points and best-fit Gaussian as line
% x=[1 2 3 4 5];y=[1 2 2.5 2 1];
% [Height,Position,Width]=gaussfit(x,y);
% plot(x,y,'o',linspace(0,8),Height.*gaussian(linspace(0,8),Position,Width))

% Copyright (c) 2012, Thomas C. O'Haver

maxy=max(y);
for p=1:length(y),
    if y(p)<(maxy/100),y(p)=maxy/100;end
end % for p=1:length(y),
z=log(y);
coef=polyfit(x,z,2);
a=coef(3);
b=coef(2);
c=coef(1);
Height=exp(a-c*(b/(2*c))^2);
Position=-b/(2*c);
Width=2.35482/(sqrt(2)*sqrt(-c));

% ----------------------------------------------------------------------
function [FitResults,GOF,baseline,coeff,BestStart,xi,yi,BootResults]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots,bipolar,minwidth,DELTA,clipheight)
% A command-line peak fitting program for time-series signals, written as a
% self-contained Matlab function in a single m-file. Uses a non-linear
% optimization algorithm to decompose a complex, overlapping-peak signal
% into its component parts. The objective is to determine whether your
% signal can be represented as the sum of fundamental underlying peaks
% shapes. Accepts signals of any length, including those with non-integer
% and non-uniform x-values. Fits any number of peaks of any of 33 curve
% shapes. This is a command line version, usable from a remote terminal. It
% is capable of making multiple trial fits with sightly different starting
% values (whose variability is controled by the 14th input argument) and
% taking the one with the lowest mean fit error (example 6). It can
% estimate the standard deviation of peak parameters from a single signal
% using the bootstrap method (example 10).
%
% Version 7: March, 2015, adds peak shapes with three unconstrained
% iterated variables: 30=voigt (variable alpha), 31=ExpGaussian (variable
% time constant), 32=Pearson (variable shape factor), 34=Gaussian/
% Lorentzian blend (variable percent). See Examples 25-28 below.
%
% Copyright (c) 2015, Thomas C. O'Haver
% 
global AA xxx PEAKHEIGHTS FIXEDPARAMETERS AUTOZERO delta BIPOLAR CLIPHEIGHT
format short g
format compact
warning off all
NumArgOut=nargout;
datasize=size(signal);
if datasize(1)<datasize(2),signal=signal';end
datasize=size(signal);
if datasize(2)==1, %  Must be isignal(Y-vector)
    X=1:length(signal); % Create an independent variable vector
    Y=signal;
else
    % Must be isignal(DataMatrix)
    X=signal(:,1); % Split matrix argument 
    Y=signal(:,2);
end
X=reshape(X,1,length(X)); % Adjust X and Y vector shape to 1 x n (rather than n x 1)
Y=reshape(Y,1,length(Y));
% If necessary, flip the data vectors so that X increases
if X(1)>X(length(X)),
    disp('X-axis flipped.')
    X=fliplr(X);
    Y=fliplr(Y);
end

% Isolate desired segment from data set for curve fitting
if nargin==1 || nargin==2,center=(max(X)-min(X))/2;window=max(X)-min(X);end
% Y=Y-min(Y);
xoffset=0;
n1=val2ind(X,center-window/2);
n2=val2ind(X,center+window/2);
if window==0,n1=1;n2=length(X);end
xx=X(n1:n2)-xoffset;
yy=Y(n1:n2);
ShapeString='Gaussian';
coeff=0;
CLIPHEIGHT=max(Y);
LOGPLOT=0;
% Define values of any missing arguments
% (signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots,bipolar,minwidth,DELTA)
switch nargin
    case 1
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        xx=X;yy=Y;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=xx(2)-xx(1);
        delta=1;
        CLIPHEIGHT=max(Y);
    case 2
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        xx=signal;yy=center;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=xx(2)-xx(1);
        delta=1;
        CLIPHEIGHT=max(Y);
    case 3
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=xx(2)-xx(1);
        delta=1;
        CLIPHEIGHT=max(Y);
    case 4 % Numpeaks specified in arguments
        peakshape=1;
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=xx(2)-xx(1);
        delta=1;
        CLIPHEIGHT=max(Y);
    case 5 % Numpeaks, peakshape specified in arguments
        extra=zeros(1,NumPeaks);
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
        CLIPHEIGHT=max(Y);
    case 6
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
    case 7
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
        CLIPHEIGHT=max(Y);
    case 8
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
        CLIPHEIGHT=max(Y);
    case 9
        AUTOZERO=autozero;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
    case 10
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
    case 11
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
        CLIPHEIGHT=max(Y);
    case 12
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        BIPOLAR=bipolar;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
        CLIPHEIGHT=max(Y);
    case 13
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        BIPOLAR=bipolar;
        MINWIDTH=minwidth;
        delta=1;
    case 14
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        BIPOLAR=bipolar;
        MINWIDTH=minwidth;
        delta=DELTA;
        CLIPHEIGHT=max(Y);
    case 15
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        BIPOLAR=bipolar;
        MINWIDTH=minwidth;
        delta=DELTA;
        CLIPHEIGHT=clipheight;
    otherwise
end % switch nargin

% Saturation Code, skips points greater than set maximum
if CLIPHEIGHT<max(Y),
    apnt=1;
    for pnt=1:length(xx),
        if yy(pnt)<CLIPHEIGHT,
            axx(apnt)=xx(pnt);
            ayy(apnt)=yy(pnt);
            apnt=apnt+1;
        end
    end
    xx=axx;yy=ayy;
end
% Default values for placeholder zeros1
if NumTrials==0;NumTrials=1;end
if isscalar(peakshape),
else
    % disp('peakshape is vector');
    shapesvector=peakshape;
    NumPeaks=length(peakshape);
    peakshape=22;
end
if peakshape==0;peakshape=1;end
if NumPeaks==0;NumPeaks=1;end
if start==0;start=calcstart(xx,NumPeaks,xoffset);end
if FIXEDPARAMETERS==0, FIXEDPARAMETERS=length(xx)/10;end
if peakshape==16;FIXEDPOSITIONS=fixedparameters;end
if peakshape==17;FIXEDPOSITIONS=fixedparameters;end
if AUTOZERO>3,AUTOZERO=3,end
if AUTOZERO<0,AUTOZERO=0,end
Heights=zeros(1,NumPeaks);
FitResults=zeros(NumPeaks,6);

% % Remove linear baseline from data segment if AUTOZERO==1
baseline=0;
bkgcoef=0;
bkgsize=round(length(xx)/10);
if bkgsize<2,bkgsize=2;end
lxx=length(xx);
if AUTOZERO==1, % linear autozero operation  
    XX1=xx(1:round(lxx/bkgsize));
    XX2=xx((lxx-round(lxx/bkgsize)):lxx);
    Y1=yy(1:(round(length(xx)/bkgsize)));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx);
    bkgcoef=polyfit([XX1,XX2],[Y1,Y2],1);  % Fit straight line to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if
if AUTOZERO==2, % Quadratic autozero operation  
    XX1=xx(1:round(lxx/bkgsize));
    XX2=xx((lxx-round(lxx/bkgsize)):lxx);
    Y1=yy(1:round(length(xx)/bkgsize));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx);
    bkgcoef=polyfit([XX1,XX2],[Y1,Y2],2);  % Fit parabola to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if autozero

PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
newstart=start;
% Assign ShapStrings
switch peakshape(1)
    case 1
        ShapeString='Gaussian';
    case 2
        ShapeString='Lorentzian';
    case 3
        ShapeString='Logistic';
    case 4
        ShapeString='Pearson';
    case 5
        ShapeString='ExpGaussian';
    case 6
        ShapeString='Equal width Gaussians';
    case 7
        ShapeString='Equal width Lorentzians';
    case 8
        ShapeString='Exp. equal width Gaussians';
    case 9
        ShapeString='Exponential Pulse';
    case 10
        ShapeString='Up Sigmoid (logistic function)';
    case 23
        ShapeString='Down Sigmoid (logistic function)';  
    case 11
        ShapeString='Fixed-width Gaussian';
    case 12
        ShapeString='Fixed-width Lorentzian';
    case 13
        ShapeString='Gaussian/Lorentzian blend';
    case 14
        ShapeString='BiGaussian';    
    case 15
        ShapeString='Breit-Wigner-Fano';   
    case 16
        ShapeString='Fixed-position Gaussians';
    case 17
        ShapeString='Fixed-position Lorentzians';
    case 18
        ShapeString='Exp. Lorentzian';
    case 19
        ShapeString='Alpha function';
    case 20
        ShapeString='Voigt (equal alphas)';
    case 21
        ShapeString='triangular';
    case 22
        ShapeString=num2str(shapesvector);
    case 24
        ShapeString='Negative Binomial Distribution';
    case 25
        ShapeString='Lognormal Distribution';
    case 26
        ShapeString='slope';
    case 27
        ShapeString='First derivative';
    case 28
        ShapeString='Polynomial';
    case 29
        ShapeString='Segmented linear';
    case 30
        ShapeString='Voigt (variable alphas)';
    case 31
        ShapeString='ExpGaussian (var. time constant)';
    case 32
        ShapeString='Pearson (var. shape constant)';
    case 33
        ShapeString='Variable Gaussian/Lorentzian';
    otherwise
end % switch peakshape
  
% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.001,'Display','off','MaxFunEvals',1000 );
LowestError=1000; % or any big number greater than largest error expected
FitParameters=zeros(1,NumPeaks.*2); 
BestStart=zeros(1,NumPeaks.*2); 
height=zeros(1,NumPeaks); 
bestmodel=zeros(size(yy));

for k=1:NumTrials, 
    % StartMatrix(k,:)=newstart;
    % disp(['Trial number ' num2str(k) ] ) % optionally prints the current trial number as progress indicator
    switch peakshape(1)
        case 1
            TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 2
            TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 3
            TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 4
            TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 5
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 6
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 7
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewlorentzian(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 8
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 9
            TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 10
            TrialParameters=fminsearch(@(lambda)(fitupsigmoid(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 23
            TrialParameters=fminsearch(@(lambda)(fitdownsigmoid(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 11
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFWGaussian(lambda,xx,yy)),fixedstart,options);
        case 12
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFWLorentzian(lambda,xx,yy)),fixedstart,options);
        case 13
            TrialParameters=fminsearch(@(lambda)(fitGL(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 14
            TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 15
            TrialParameters=fminsearch(@(lambda)(fitBWF(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 16
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
                fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc);
            end
            TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(Peak)<MINWIDTH,
                    TrialParameters(Peak)=MINWIDTH;
                end
            end
        case 17
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
                fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc);
            end
            TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(Peak)<MINWIDTH,
                    TrialParameters(Peak)=MINWIDTH;
                end
            end
        case 18
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexplorentzian(lambda,zxx,zyy,-extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 19
            TrialParameters=fminsearch(@(lambda)(fitalphafunction(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 20
            TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 21
            TrialParameters=fminsearch(@(lambda)(fittriangular(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 22
            TrialParameters=fminsearch(@(lambda)(fitmultiple(lambda,xx,yy,shapesvector,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH(Peak),
                    TrialParameters(2*Peak)=MINWIDTH(Peak);
                end
            end
        case 24
            TrialParameters=fminsearch(@(lambda)(fitnbinpdf(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 25
            TrialParameters=fminsearch(@(lambda)(fitlognpdf(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 26
            TrialParameters=fminsearch(@(lambda)(fitlinslope(lambda,xx,yy)),polyfit(xx,yy,1),options);
             coeff=TrialParameters;
        case 27
            TrialParameters=fminsearch(@(lambda)(fitd1gauss(lambda,xx,yy)),newstart,options);
        case 28
            coeff=fitpolynomial(xx,yy,extra);
            TrialParameters=coeff;
        case 29
            cnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cnewstart(pc)=newstart(2.*pc-1)+(delta*(rand-.5)/50);
            end
            TrialParameters=fminsearch(@(lambda)(fitsegmented(lambda,xx,yy)),cnewstart,options);
        case 30
            nn=max(xx)-min(xx);
            start=[];
            startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker)+ xoffset;
                start=[start markx nn/5 extra];
            end % for marker
             newstart=start;
            for parameter=1:3:3*NumPeaks,
                newstart(parameter)=newstart(parameter)*(1+randn/100);
                newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
            end
            TrialParameters=fminsearch(@(lambda)(fitvoigtv(lambda,xx,yy)),newstart);
         case 31
            nn=max(xx)-min(xx);
            start=[];
            startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker)+ xoffset;
                start=[start markx nn/5 extra];
            end % for marker
             newstart=start;
            for parameter=1:3:3*NumPeaks,
                newstart(parameter)=newstart(parameter)*(1+randn/100);
                newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
            end
            % newstart=newstart
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexpgaussianv(lambda,zxx,zyy)),newstart);
        case 32
            nn=max(xx)-min(xx);
            start=[];
            startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker)+ xoffset;
                start=[start markx nn/5 extra];
            end % for marker
             newstart=start;
            for parameter=1:3:3*NumPeaks,
                newstart(parameter)=newstart(parameter)*(1+randn/100);
                newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
            end
            % newstart=newstart
            TrialParameters=fminsearch(@(lambda)(fitpearsonv(lambda,xx,yy)),newstart);
        case 33
             nn=max(xx)-min(xx);
            start=[];
            startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker)+ xoffset;
                start=[start markx nn/5 extra];
            end % for marker
             newstart=start;
            for parameter=1:3:3*NumPeaks,
                newstart(parameter)=newstart(parameter)*(1+randn/100);
                newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
            end
            % newstart=newstart
            TrialParameters=fminsearch(@(lambda)(fitGLv(lambda,xx,yy)),newstart);
        otherwise
    end % switch peakshape

% Construct model from Trial parameters
A=zeros(NumPeaks,n);
for m=1:NumPeaks,
    switch peakshape(1)
        case 1
            A(m,:)=gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 2
            A(m,:)=lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 3
            A(m,:)=logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 4
            A(m,:)=pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 5
            A(m,:)=expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
        case 6
            A(m,:)=gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
        case 7
            A(m,:)=lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
        case 8
            A(m,:)=expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)';
        case 9
            A(m,:)=exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 10
            A(m,:)=upsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 11
            A(m,:)=gaussian(xx,TrialParameters(m),FIXEDPARAMETERS);
        case 12
            A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDPARAMETERS);
        case 13
            A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 14
            A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 15
            A(m,:)=BWF(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
        case 16
            A(m,:)=gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
        case 17
            A(m,:)=lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
        case 18
            A(m,:)=explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
        case 19
            A(m,:)=alphafunction(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 20
            A(m,:)=voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
        case 21
            A(m,:)=triangular(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 22
            A(m,:)=peakfunction(shapesvector(m),xx,TrialParameters(2*m-1),TrialParameters(2*m),extra(m));        
        case 23
            A(m,:)=downsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));        
        case 24
            A(m,:)=nbinpdf(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 25
            A(m,:)=lognormal(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 26
            A(m,:)=linslope(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 27
            A(m,:)=d1gauss(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 28
            A(m,:)=polynomial(xx,coeff);
        case 29
            A(m,:)=segmented(xx,yy,PEAKHEIGHTS);
        case 30
            A(m,:)=voigt(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
        case 31
            A(m,:)=expgaussian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),-TrialParameters(3*m));        
        case 32
            A(m,:)=pearson(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
        case 33
            A(m,:)=GL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
        otherwise
    end % switch
    for parameter=1:2:2*NumPeaks,
        newstart(parameter)=newstart(parameter)*(1+delta*(rand-.5)/500);
        newstart(parameter+1)=newstart(parameter+1)*(1+delta*(rand-.5)/100);
    end
end % for NumPeaks

% Multiplies each row by the corresponding amplitude and adds them up
if peakshape(1)==29, % Segmented linear
    model=segmented(xx,yy,PEAKHEIGHTS);
    TrialParameters=PEAKHEIGHTS;
    Heights=ones(size(PEAKHEIGHTS));
else
    if AUTOZERO==3,
        baseline=PEAKHEIGHTS(1);
        Heights=PEAKHEIGHTS(2:1+NumPeaks);
        model=Heights'*A+baseline;
    else
%          size(PEAKHEIGHTS) % error check
%          size(A)
        model=PEAKHEIGHTS'*A;
        Heights=PEAKHEIGHTS;
        baseline=0;
    end
end
if peakshape(1)==28, % polynomial;
    model=polynomial(xx,coeff);
    TrialParameters=PEAKHEIGHTS;
    Heights=ones(size(PEAKHEIGHTS));
end
% Compare trial model to data segment and compute the fit error
    MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
  % Take only the single fit that has the lowest MeanFitError
  if MeanFitError<LowestError, 
      if min(Heights)>=-BIPOLAR*10^100,  % Consider only fits with positive peak heights
        LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
        FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
        BestStart=newstart; % Assign BestStart to the start with the lowest MeanFitError
        height=Heights; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
        bestmodel=model; % Assign bestmodel to the model with the lowest MeanFitError
      end % if min(PEAKHEIGHTS)>0
  end % if MeanFitError<LowestError
%  ErrorVector(k)=MeanFitError;
end % for k (NumTrials)
    Rsquared=1-(norm(yy-bestmodel)./norm(yy-mean(yy)));
    SStot=sum((yy-mean(yy)).^2);
    SSres=sum((yy-bestmodel).^2);
    Rsquared=1-(SSres./SStot);
    GOF=[LowestError Rsquared];
% Uncomment following 4 lines to monitor trail fit starts and errors.
% StartMatrix=StartMatrix;
% ErrorVector=ErrorVector;
% matrix=[StartMatrix ErrorVector']
% std(StartMatrix)
% Construct model from best-fit parameters
AA=zeros(NumPeaks,600);
xxx=linspace(min(xx),max(xx),600);
% xxx=linspace(min(xx)-length(xx),max(xx)+length(xx),200);
for m=1:NumPeaks,
   switch peakshape(1)
    case 1
        AA(m,:)=gaussian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 2
        AA(m,:)=lorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 3
        AA(m,:)=logistic(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 4
        AA(m,:)=pearson(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 5
        AA(m,:)=expgaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))';
    case 6
        AA(m,:)=gaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
    case 7
        AA(m,:)=lorentzian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
    case 8
        AA(m,:)=expgaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1),-extra*length(xxx)./length(xx))';
    case 9
        AA(m,:)=exppulse(xxx,FitParameters(2*m-1),FitParameters(2*m));  
    case 10
        AA(m,:)=upsigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m));   
    case 11
        AA(m,:)=gaussian(xxx,FitParameters(m),FIXEDPARAMETERS);
    case 12
        AA(m,:)=lorentzian(xxx,FitParameters(m),FIXEDPARAMETERS);
    case 13
        AA(m,:)=GL(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 14
        AA(m,:)=BiGaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 15
        AA(m,:)=BWF(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 16
        AA(m,:)=gaussian(xxx,FIXEDPOSITIONS(m),FitParameters(m));
    case 17
        AA(m,:)=lorentzian(xxx,FIXEDPOSITIONS(m),FitParameters(m));
    case 18
        AA(m,:)=explorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))';
    case 19
        AA(m,:)=alphafunction(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 20
        AA(m,:)=voigt(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 21
        AA(m,:)=triangular(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 22
        AA(m,:)=peakfunction(shapesvector(m),xxx,FitParameters(2*m-1),FitParameters(2*m),extra(m));        
    case 23
        AA(m,:)=downsigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m));  
    case 24
        AA(m,:)=nbinpdf(xxx,FitParameters(2*m-1),FitParameters(2*m));    
    case 25
        AA(m,:)=lognormal(xxx,FitParameters(2*m-1),FitParameters(2*m));    
    case 26
        AA(m,:)=linslope(xxx,FitParameters(2*m-1),FitParameters(2*m));   
    case 27
        AA(m,:)=d1gauss(xxx,FitParameters(2*m-1),FitParameters(2*m));  
    case 28
        AA(m,:)=polynomial(xxx,coeff);
    case 29
    case 30
        AA(m,:)=voigt(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m));        
    case 31
        AA(m,:)=expgaussian(xxx,FitParameters(3*m-2),FitParameters(3*m-1),-FitParameters(3*m)*length(xxx)./length(xx));        
    case 32
        AA(m,:)=pearson(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m));        
    case 33
        AA(m,:)=GL(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m)); 
       otherwise
  end % switch
end % for NumPeaks

% Multiplies each row by the corresponding amplitude and adds them up
if peakshape(1)==29, % Segmented linear
    mmodel=segmented(xx,yy,PEAKHEIGHTS);
    baseline=0;
else
    heightsize=size(height');
    AAsize=size(AA);
    if heightsize(2)==AAsize(1),
        mmodel=height'*AA+baseline;
    else
        mmodel=height*AA+baseline;
    end
end
% Top half of the figure shows original signal and the fitted model.
if plots,
    subplot(2,1,1);plot(xx+xoffset,yy,'b.'); % Plot the original signal in blue dots
    hold on
end
if peakshape(1)==28, % Polynomial
     yi=polynomial(xxx,coeff);
else
    for m=1:NumPeaks,
        if plots, plot(xxx+xoffset,height(m)*AA(m,:)+baseline,'g'),end  % Plot the individual component peaks in green lines
        area(m)=trapz(xxx+xoffset,height(m)*AA(m,:)); % Compute the area of each component peak using trapezoidal method
        yi(m,:)=height(m)*AA(m,:); % Place y values of individual model peaks into matrix yi
    end
end
xi=xxx+xoffset; % Place the x-values of the individual model peaks into xi

if plots,
    % Mark starting peak positions with vertical dashed magenta lines
    if peakshape(1)==16||peakshape(1)==17
    else
        if peakshape(1)==29, % Segmented linear
            subplot(2,1,1);plot([PEAKHEIGHTS' PEAKHEIGHTS'],[0 max(yy)],'m--')
        else
            for marker=1:NumPeaks,
                markx=BestStart((2*marker)-1);
                subplot(2,1,1);plot([markx+xoffset markx+xoffset],[0 max(yy)],'m--')
            end % for
        end
    end % if peakshape

    % Plot the total model (sum of component peaks) in red lines
    if peakshape(1)==29, % Segmented linear
        mmodel=segmented(xx,yy,PEAKHEIGHTS);
       plot(xx+xoffset,mmodel,'r');  
    else
       plot(xxx+xoffset,mmodel,'r');  
    end
    hold off;
    lyy=min(yy);
    uyy=max(yy)+(max(yy)-min(yy))/10;
    if BIPOLAR,
        axis([min(xx) max(xx) lyy uyy]);
        ylabel('+ - mode')
    else
        axis([min(xx) max(xx) 0 uyy]);
        ylabel('+ mode')
    end
    switch AUTOZERO,
        case 0
            title(['peakfit.m Version 7   No baseline correction'])
        case 1
            title(['peakfit.m Version 7   Linear baseline subtraction'])
        case 2
            title(['peakfit.m Version 7   Quadratic subtraction baseline'])
        case 3
            title(['peakfit.m Version 7   Flat baseline correction'])
    end
 
    switch peakshape(1)
        case {4,20}
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      Shape Constant = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000) ] )
        case {5,8,18}
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      Time Constant = ' num2str(extra)   '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000)  ] )
        case 13
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      % Gaussian = ' num2str(extra)   '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000)  ] )
        case {14,15,22}
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      extra = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000) ] )
        case 28
            xlabel(['Shape = ' ShapeString '      Order = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '%  R2 = ' num2str(round(1000*LowestError)/1000) ] )
        otherwise
            if peakshape(1)==29, % Segmented linear
                xlabel(['Breakpoints = ' num2str(NumPeaks) '     Shape = ' ShapeString  '     Error = ' num2str(round(1000*LowestError)/1000) '%  R2 = ' num2str(round(100000*Rsquared)/100000) ] )
            else
                xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH)  '     Error = ' num2str(round(1000*LowestError)/1000) '%  R2 = ' num2str(round(100000*Rsquared)/100000) ] )
            end % if peakshape(1)==29
    end % switch peakshape(1)

    % Bottom half of the figure shows the residuals and displays RMS error
    % between original signal and model
    residual=yy-bestmodel;
    subplot(2,1,2);plot(xx+xoffset,residual,'r.')
    axis([min(xx)+xoffset max(xx)+xoffset min(residual) max(residual)]);
    xlabel('Residual Plot')
    if NumTrials>1,
       title(['Best of ' num2str(NumTrials) ' fits'])
    else
       title(['Single fit'])
    end
end % if plots

% Put results into a matrix FitResults, one row for each peak, showing peak index number,
% position, amplitude, and width.
FitResults=zeros(NumPeaks,6);
switch peakshape(1),
    case {6,7,8}, % equal-width peak models only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
            end
        end
    case {11,12}, % Fixed-width shapes only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS area(m)]];
            end
        end
    case {16,17}, % Fixed-position shapes only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)];
            else
                FitResults=[FitResults ; [round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]];
            end
        end
    case 28,   % Simple polynomial fit
        FitResults=PEAKHEIGHTS;
    case 29, % Segmented linear fit
        FitResults=PEAKHEIGHTS;
    case {30,31,32,33} % Special case of shapes with 3 iterated variables
        for m=1:NumPeaks,
            if m==1,
                FitResults=[round(m) FitParameters(3*m-2) height(m) abs(FitParameters(3*m-1)) area(m) FitParameters(3*m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(3*m-2) height(m) abs(FitParameters(3*m-1)) area(m)] FitParameters(3*m)];
            end
        end
    otherwise % Normal shapes with 2 iterated variables
        for m=1:NumPeaks,
            if m==1,
                FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]];
            end % if m==1
        end % for m=1:NumPeaks,
end % switch peakshape(1)
  
% Display Fit Results on lower graph
if plots,
    % Display Fit Results on lower  graph
    subplot(2,1,2);
    startx=min(xx)+(max(xx)-min(xx))./20;
    dxx=(max(xx)-min(xx))./10;
    dyy=((max(residual)-min(residual))./10);
    starty=max(residual)-dyy;
    FigureSize=get(gcf,'Position');
    switch peakshape(1)
        case {9,19,11}  % Pulse and sigmoid shapes only
            text(startx,starty+dyy/2,['Peak #          tau1           Height           tau2             Area'] );
        case 28, % Polynomial
            text(startx,starty+dyy/2,['Polynomial coefficients'] );
        case 29 % Segmented linear
             text(startx,starty+dyy/2,['x-axis breakpoints'] );
        case {30,31,32,33} % Special case of shapes with 3 iterated variables
            text(startx,starty+dyy/2,['Peak #          Position        Height         Width             Area       Shape factor'] );            
        otherwise
            text(startx,starty+dyy/2,['Peak #          Position        Height         Width             Area '] );
    end
    % Display FitResults using sprintf
    if peakshape(1)==28||peakshape(1)==29, % Polynomial or segmented linear
        for number=1:length(FitResults),
            column=1;
            itemstring=sprintf('%0.4g',FitResults(number));
            xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)));
            yposition=starty-number.*dyy.*(400./FigureSize(4));
            text(xposition,yposition,['                ' itemstring]);
        end
    else
        for peaknumber=1:NumPeaks,
            for column=1:5,
                itemstring=sprintf('%0.4g',FitResults(peaknumber,column));
                xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)));
                yposition=starty-peaknumber.*dyy.*(400./FigureSize(4));
                text(xposition,yposition,itemstring);
            end
        end
        xposition=startx;
        yposition=starty-(peaknumber+1).*dyy.*(400./FigureSize(4));
        if AUTOZERO==3,
            text(xposition,yposition,[ 'Baseline= ' num2str(baseline) ]);
        end % if AUTOZERO
    end % if peakshape(1)
    if peakshape(1)==30 || peakshape(1)==31 || peakshape(1)==32 || peakshape(1)==33,
        for peaknumber=1:NumPeaks,
            column=6;
            itemstring=sprintf('%0.4g',FitParameters(3*peaknumber));
            xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)));
            yposition=starty-peaknumber.*dyy.*(400./FigureSize(4));
            text(xposition,yposition,itemstring);
        end
    end
end % if plots

if NumArgOut==8,
    if plots,disp('Computing bootstrap sampling statistics.....'),end
    BootstrapResultsMatrix=zeros(6,100,NumPeaks);
    BootstrapErrorMatrix=zeros(1,100,NumPeaks);
    clear bx by
    tic;
    for trial=1:100,
        n=1;
        bx=xx;
        by=yy;
        while n<length(xx)-1,
            if rand>.5,
                bx(n)=xx(n+1);
                by(n)=yy(n+1);
            end
            n=n+1;
        end
        bx=bx+xoffset;
        [FitResults,BootFitError]=fitpeaks(bx,by,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,FIXEDPARAMETERS);
        for peak=1:NumPeaks,
            switch peakshape(1)
                case {30,31,32,33}
                    BootstrapResultsMatrix(1:6,trial,peak)=FitResults(peak,1:6);
                otherwise
                    BootstrapResultsMatrix(1:5,trial,peak)=FitResults(peak,1:5);
            end
            BootstrapErrorMatrix(:,trial,peak)=BootFitError;
        end
    end
    if plots,toc;end
    for peak=1:NumPeaks,
        if plots,
            disp(' ')
            disp(['Peak #',num2str(peak) '         Position    Height       Width       Area      Shape Factor']);
        end % if plots
        BootstrapMean=mean(real(BootstrapResultsMatrix(:,:,peak)'));
        BootstrapSTD=std(BootstrapResultsMatrix(:,:,peak)');
        BootstrapIQR=iqr(BootstrapResultsMatrix(:,:,peak)');
        PercentRSD=100.*BootstrapSTD./BootstrapMean;
        PercentIQR=100.*BootstrapIQR./BootstrapMean;
        BootstrapMean=BootstrapMean(2:6);
        BootstrapSTD=BootstrapSTD(2:6);
        BootstrapIQR=BootstrapIQR(2:6);
        PercentRSD=PercentRSD(2:6);
        PercentIQR=PercentIQR(2:6);
        if plots,
            disp(['Bootstrap Mean: ', num2str(BootstrapMean)])
            disp(['Bootstrap STD:  ', num2str(BootstrapSTD)])
            disp(['Bootstrap IQR:  ', num2str(BootstrapIQR)])
            disp(['Percent RSD:    ', num2str(PercentRSD)])
            disp(['Percent IQR:    ', num2str(PercentIQR)])
        end % if plots
        BootResults(peak,:)=[BootstrapMean BootstrapSTD PercentRSD BootstrapIQR PercentIQR];
    end % peak=1:NumPeaks,
end % if NumArgOut==8,
if AUTOZERO==3;
else
    baseline=bkgcoef;
end
% ----------------------------------------------------------------------
function [FitResults,LowestError]=fitpeaks(xx,yy,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,fixedparameters)
% Based on peakfit Version 3: June, 2012. 
global PEAKHEIGHTS FIXEDPARAMETERS AUTOZERO BIPOLAR MINWIDTH coeff
format short g
format compact
warning off all
FIXEDPARAMETERS=fixedparameters;
xoffset=0;
if start==0;start=calcstart(xx,NumPeaks,xoffset);end
PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
newstart=start;
coeff=0;
LOGPLOT=0;

% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.001,'Display','off','MaxFunEvals',1000 );
LowestError=1000; % or any big number greater than largest error expected
FitParameters=zeros(1,NumPeaks.*2); 
BestStart=zeros(1,NumPeaks.*2); 
height=zeros(1,NumPeaks); 
bestmodel=zeros(size(yy));
for k=1:NumTrials,
    % StartVector=newstart
    switch peakshape(1)
        case 1
            TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 2
            TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 3
            TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 4
            TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 5
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 6
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 7
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewlorentzian(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 8
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 9
            TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 10
            TrialParameters=fminsearch(@(lambda)(fitupsigmoid(lambda,xx,yy)),newstar,optionst);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 11
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFWGaussian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 12
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFWLorentzian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 13
            TrialParameters=fminsearch(@(lambda)(fitGL(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 14
            TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 15
            TrialParameters=fminsearch(@(lambda)(fitBWF(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 16
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(Peak)<MINWIDTH,
                    TrialParameters(Peak)=MINWIDTH;
                end
            end
        case 17
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(Peak)<MINWIDTH,
                    TrialParameters(Peak)=MINWIDTH;
                end
            end
        case 18
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexplorentzian(lambda,zxx,zyy,-extra)),newstart,options);
        case 19
            TrialParameters=fminsearch(@(lambda)(alphafunction(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 20
            TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 21
            TrialParameters=fminsearch(@(lambda)(fittriangular(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 22
            TrialParameters=fminsearch(@(lambda)(fitmultiple(lambda,xx,yy,shapesvector,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH(Peak),
                    TrialParameters(2*Peak)=MINWIDTH(Peak);
                end
            end
        case 23
            TrialParameters=fminsearch(@(lambda)(fitdownsigmoid(lambda,xx,yy)),newstar,optionst);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 24
            TrialParameters=fminsearch(@(lambda)(fitnbinpdf(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 25
            TrialParameters=fminsearch(@(lambda)(fitlognpdf(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 26
            TrialParameters=fminsearch(@(lambda)(fitlinslope(lambda,xx,yy)),polyfit(xx,yy,1),options);
        coeff=TrialParameters;
        case 27
            TrialParameters=fminsearch(@(lambda)(fitd1gauss(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 28
            TrialParameters=fitpolynomial(xx,yy,extra);
        case 29
            TrialParameters=fminsearch(@(lambda)(fitsegmented(lambda,xx,yy)),newstart,options);
        case 30
            TrialParameters=fminsearch(@(lambda)(fitvoigtv(lambda,xx,yy)),newstart);
        case 31
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexpgaussianv(lambda,zxx,zyy)),newstart);
        case 32
            TrialParameters=fminsearch(@(lambda)(fitpearsonv(lambda,xx,yy)),newstart);
        case 33
            TrialParameters=fminsearch(@(lambda)(fitGLv(lambda,xx,yy)),newstart);
        otherwise
    end % switch peakshape
    
for peaks=1:NumPeaks,
     peakindex=2*peaks-1;
     newstart(peakindex)=start(peakindex)-xoffset;
end

    % Construct model from Trial parameters
    A=zeros(NumPeaks,n);
    for m=1:NumPeaks,
        switch peakshape(1)
            case 1
                A(m,:)=gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 2
                A(m,:)=lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 3
                A(m,:)=logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 4
                A(m,:)=pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 5
                A(m,:)=expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
            case 6
                A(m,:)=gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
            case 7
                A(m,:)=lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
            case 8
                A(m,:)=expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)';
            case 9
                A(m,:)=exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 10
                A(m,:)=upsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 11
                A(m,:)=gaussian(xx,TrialParameters(m),FIXEDPARAMETERS);
            case 12
                A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDPARAMETERS);
            case 13
                A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 14
                A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 15
                A(m,:)=BWF(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 16
                A(m,:)=gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
            case 17
                A(m,:)=lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
            case 18
                A(m,:)=explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
            case 19
                A(m,:)=alphafunction(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 20
                A(m,:)=voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 21
                A(m,:)=triangular(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 22
                A(m,:)=peakfunction(shapesvector(m),xx,TrialParameters(2*m-1),TrialParameters(2*m),extra(m));
            case 23
                A(m,:)=downsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));      
            case 24
                A(m,:)=nbinpdf(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 25
                A(m,:)=lognormal(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 26
                A(m,:)=linslope(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 27
                A(m,:)=d1gauss(xx,TrialParameters(2*m-1),TrialParameters(2*m));       
            case 28
                A(m,:)=polynomial(xx,TrialParameters(2*m-1),TrialParameters(2*m));       
            case 29
                A(m,:)=segmented(xx,yy,PEAKHEIGHTS);
            case 30
                A(m,:)=voigt(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
            case 31
                A(m,:)=expgaussian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
            case 32
                A(m,:)=pearson(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
            case 33
                A(m,:)=GL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));
        end % switch
    end % for
    
    % Multiplies each row by the corresponding amplitude and adds them up
    if peakshape(1)==29, % Segmented linear
        model=segmented(xx,yy,PEAKHEIGHTS);
        TrialParameters=coeff;
        Heights=ones(size(coeff));
    else
        if AUTOZERO==3,
            baseline=PEAKHEIGHTS(1);
            Heights=PEAKHEIGHTS(2:1+NumPeaks);
            model=Heights'*A+baseline;
        else
            model=PEAKHEIGHTS'*A;
            Heights=PEAKHEIGHTS;
            baseline=0;
        end
    end
    
    % Compare trial model to data segment and compute the fit error
    MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
    % Take only the single fit that has the lowest MeanFitError
    if MeanFitError<LowestError,
        if min(Heights)>=-BIPOLAR*10^100,  % Consider only fits with positive peak heights
            LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
            FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
            height=Heights; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
        end % if min(PEAKHEIGHTS)>0
    end % if MeanFitError<LowestError
end % for k (NumTrials)
    Rsquared=1-(norm(yy-bestmodel)./norm(yy-mean(yy)));
    SStot=sum((yy-mean(yy)).^2);
    SSres=sum((yy-bestmodel).^2);
    Rsquared=1-(SSres./SStot);
    GOF=[LowestError Rsquared];
for m=1:NumPeaks,
    area(m)=trapz(xx+xoffset,height(m)*A(m,:)); % Compute the area of each component peak using trapezoidal method
end

% Put results into a matrix FitResults, one row for each peak, showing peak index number,
% position, amplitude, and width.
FitResults=zeros(NumPeaks,6);
switch peakshape(1),
    case {6,7,8}, % equal-width peak models only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
            end
        end
    case {11,12}, % Fixed-width shapes only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS area(m)]];
            end
        end
    case {16,17}, % Fixed-position shapes only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)];
            else
                FitResults=[FitResults ; [round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]];
            end
        end
    case 28,   % Simple polynomial fit
        FitResults=PEAKHEIGHTS;
    case 29, % Segmented linear fit
        FitResults=PEAKHEIGHTS;
    case {30,31,32,33} % Special case of shapes with 3 iterated variables
        for m=1:NumPeaks,
            if m==1,
                FitResults=[round(m) FitParameters(3*m-2) height(m) abs(FitParameters(3*m-1)) area(m) FitParameters(3*m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(3*m-2) height(m) abs(FitParameters(3*m-1)) area(m) FitParameters(3*m)]];
            end
        end
    otherwise % Normal shapes with 2 iterated variables
        for m=1:NumPeaks,
            if m==1,
                FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]];
            end % if m==1
        end % for m=1:NumPeaks,
end % switch peakshape(1)
% ----------------------------------------------------------------------
function start=calcstart(xx,NumPeaks,xoffset)
  n=max(xx)-min(xx);
  start=[];
  startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
  for marker=1:NumPeaks,
      markx=startpos(marker)+ xoffset;
      start=[start markx n/ (3.*NumPeaks)];
  end % for marker
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);
% ----------------------------------------------------------------------
function err = fitgaussian(lambda,t,y)
% Fitting function for a Gaussian band signal.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
numpeaks=round(length(lambda)/2);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
%    if lambda(2*j)<MINWIDTH,lambda(2*j)=MINWIDTH;end
    A(:,j) = gaussian(t,lambda(2*j-1),lambda(2*j))';
end 
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitewgaussian(lambda,t,y)
% Fitting function for a Gaussian band signal with equal peak widths.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),lambda(numpeaks+1))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = FitFWGaussian(lambda,t,y)
%	Fitting function for a fixed width Gaussian
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),FIXEDPARAMETERS)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = FitFPGaussian(lambda,t,y)
%	Fitting function for fixed-position Gaussians
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,FIXEDPARAMETERS(j), lambda(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = FitFPLorentzian(lambda,t,y)
%	Fitting function for fixed-position Lorentzians
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,FIXEDPARAMETERS(j), lambda(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFWLorentzian(lambda,t,y)
%	Fitting function for fixed width Lorentzian
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),FIXEDPARAMETERS)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitewlorentzian(lambda,t,y)
% Fitting function for a Lorentzian band signal with equal peak widths.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),lambda(numpeaks+1))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);
% ----------------------------------------------------------------------
function err = fitlorentzian(lambda,t,y)
%	Fitting function for single lorentzian, lambda(1)=position, lambda(2)=width
%	Fitgauss assumes a lorentzian function 
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = lorentzian(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = lorentzian(x,position,width)
% lorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix
% position and width scalar
% T. C. O'Haver, 1988
% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);
% ----------------------------------------------------------------------
function err = fitlogistic(lambda,t,y)
%	Fitting function for logistic, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fitlogistic assumes a logistic function 
%  T. C. O'Haver, May 2006
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = logistic(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = logistic(x,pos,wid)
% logistic function.  pos=position; wid=half-width (both scalar)
% logistic(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991 
n = exp(-((x-pos)/(.477.*wid)) .^2);
g = (2.*n)./(1+n);
% ----------------------------------------------------------------------
function err = fittriangular(lambda,t,y)
%	Fitting function for triangular, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fittriangular assumes a triangular function 
%  T. C. O'Haver, May 2006
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = triangular(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = triangular(x,pos,wid)
%triangle function.  pos=position; wid=half-width (both scalar)
%trianglar(x,pos,wid), where x may be scalar or vector,
%pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991
% Example
% x=[0:.1:10];plot(x,trianglar(x,5.5,2.3),'.')
g=1-(1./wid) .*abs(x-pos);
for i=1:length(x),  
if g(i)<0,g(i)=0;end
end
% ----------------------------------------------------------------------
function err = fitpearson(lambda,t,y,shapeconstant)
%   Fitting functions for a Pearson 7 band signal.
% T. C. O'Haver (toh@umd.edu),   Version 1.3, October 23, 2006.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = pearson(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitpearsonv(lambda,t,y)
% Fitting functions for pearson function with independently variable
% percent Gaussian
% T. C. O'Haver (toh@umd.edu), 2015.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = pearson(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = pearson(x,pos,wid,m)
% Pearson VII function. 
% g = pearson7(x,pos,wid,m) where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% m=some number
%  T. C. O'Haver, 1990  
g=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m;
% ----------------------------------------------------------------------
function err = fitexpgaussian(lambda,t,y,timeconstant)
%   Fitting functions for a exponentially-broadened Gaussian band signal.
%  T. C. O'Haver, October 23, 2006.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = expgaussian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitexplorentzian(lambda,t,y,timeconstant)
%   Fitting functions for a exponentially-broadened lorentzian band signal.
%  T. C. O'Haver, 2013.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = explorentzian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitexpewgaussian(lambda,t,y,timeconstant)
% Fitting function for exponentially-broadened Gaussian bands with equal peak widths.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = expgaussian(t,lambda(j),lambda(numpeaks+1),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitexpgaussianv(lambda,t,y)
% Fitting functions for  exponentially-broadened Gaussians with
% independently variable time constants
% T. C. O'Haver (toh@umd.edu), 2015.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = expgaussian(t,lambda(3*j-2),lambda(3*j-1),-lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = expgaussian(x,pos,wid,timeconstant)
%  Exponentially-broadened gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 2006
g = exp(-((x-pos)./(0.6005615.*wid)) .^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function g = explorentzian(x,pos,wid,timeconstant)
%  Exponentially-broadened lorentzian(x,pos,wid) = lorentzian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 2013
g = ones(size(x))./(1+((x-pos)./(0.5.*wid)).^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function yb = ExpBroaden(y,t)
% ExpBroaden(y,t) zero pads y and convolutes result by an exponential decay
% of time constant t by multiplying Fourier transforms and inverse
% transforming the result.
hly=round(length(y)./2);
ey=[y(1).*ones(1,hly)';y;y(length(y)).*ones(1,hly)'];
% figure(2);plot(ey);figure(1);
fy=fft(ey);
a=exp(-(1:length(fy))./t);
fa=fft(a);
fy1=fy.*fa';
ybz=real(ifft(fy1))./sum(a);
yb=ybz(hly+2:length(ybz)-hly+1);
% ----------------------------------------------------------------------
function err = fitexppulse(tau,x,y)
% Iterative fit of the sum of exponential pulses
% of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = exppulse(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = exppulse(x,t1,t2)
% Exponential pulse of the form 
% g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
e=(x-t1)./t2;
p = 4*exp(-e).*(1-exp(-e));
p=p .* (p>0);
g = p';
% ----------------------------------------------------------------------
function err = fitalphafunction(tau,x,y)
% Iterative fit of the sum of alpha funciton
% of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = alphafunction(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = alphafunction(x,pos,spoint)
% alpha function.  pos=position; wid=half-width (both scalar)
% alphafunction(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% Taekyung Kwon, July 2013  
g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
for m=1:length(x);if g(m)<0;g(m)=0;end;end
% ----------------------------------------------------------------------
function err = fitdownsigmoid(tau,x,y)
% Fitting function for iterative fit to the sum of
% downward moving sigmiods 
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = downsigmoid(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitupsigmoid(tau,x,y)
% Fitting function for iterative fit to the sum of
% upwards moving sigmiods
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = upsigmoid(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g=downsigmoid(x,t1,t2)
 % down step sigmoid
g=.5-.5*erf(real((x-t1)/sqrt(2*t2)));
% ----------------------------------------------------------------------
function g=upsigmoid(x,t1,t2)
% up step sigmoid
g=1/2 + 1/2* erf(real((x-t1)/sqrt(2*t2))); 
% ----------------------------------------------------------------------
function err = fitGL(lambda,t,y,shapeconstant)
%   Fitting functions for Gaussian/Lorentzian blend.
% T. C. O'Haver (toh@umd.edu), 2012.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = GL(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitGLv(lambda,t,y)
% Fitting functions for Gaussian/Lorentzian blend function with
% independently variable percent Gaussian
% T. C. O'Haver (toh@umd.edu), 2015.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = GL(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = GL(x,pos,wid,m)
% Gaussian/Lorentzian blend. m = percent Gaussian character
% pos=position; wid=half-width
% m = percent Gaussian character.
%  T. C. O'Haver, 2012
% sizex=size(x)
% sizepos=size(pos)
% sizewid=size(wid)
% sizem=size(m)
g=2.*((m/100).*gaussian(x,pos,wid)+(1-(m(1)/100)).*lorentzian(x,pos,wid))/2;
% ----------------------------------------------------------------------
function err = fitvoigt(lambda,t,y,shapeconstant)
% Fitting functions for Voigt profile function
% T. C. O'Haver (toh@umd.edu), 2013.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = voigt(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitvoigtv(lambda,t,y)
% Fitting functions for Voigt profile function with independently variable
% alphas
% T. C. O'Haver (toh@umd.edu), 2015.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = voigt(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g=voigt(xx,pos,gD,alpha)
% Voigt profile function. xx is the independent variable (energy,
% wavelength, etc), gD is the Doppler (Gaussian) width, and alpha is the
% shape constant (ratio of the Lorentzian width gL to the Doppler width gD.
% Based on Chong Tao's "Voigt lineshape spectrum simulation", 
% File ID: #26707
% alpha=alpha
gL=alpha.*gD;
gV = 0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2);
x = gL/gV;
y = abs(xx-pos)/gV;
g = 1/(2*gV*(1.065 + 0.447*x + 0.058*x^2))*((1-x)*exp(-0.693.*y.^2) + (x./(1+y.^2)) + 0.016*(1-x)*x*(exp(-0.0841.*y.^2.25)-1./(1 + 0.021.*y.^2.25)));
g=g./max(g);
% ----------------------------------------------------------------------
function err = fitBiGaussian(lambda,t,y,shapeconstant)
%   Fitting functions for BiGaussian.
% T. C. O'Haver (toh@umd.edu),  2012.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = BiGaussian(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = BiGaussian(x,pos,wid,m)
% BiGaussian (different widths on leading edge and trailing edge).
% pos=position; wid=width 
% m determines shape; symmetrical if m=50.
%  T. C. O'Haver, 2012
lx=length(x);
hx=val2ind(x,pos);
g(1:hx)=gaussian(x(1:hx),pos,wid*(m/100));
g(hx+1:lx)=gaussian(x(hx+1:lx),pos,wid*(1-m/100));
% ----------------------------------------------------------------------
function err = fitBWF(lambda,t,y,shapeconstant)
%   Fitting function for Breit-Wigner-Fano.
% T. C. O'Haver (toh@umd.edu),  2014.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = BWF(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = BWF(x,pos,wid,m)
% BWF (Breit-Wigner-Fano) http://en.wikipedia.org/wiki/Fano_resonance
% pos=position; wid=width; m=Fano factor
%  T. C. O'Haver, 2014
y=((m*wid/2+x-pos).^2)./(((wid/2).^2)+(x-pos).^2);
% y=((1+(x-pos./(m.*wid))).^2)./(1+((x-pos)./wid).^2);
g=y./max(y);
% ----------------------------------------------------------------------
function err = fitnbinpdf(tau,x,y)
% Fitting function for iterative fit to the sum of
% Negative Binomial Distributions
% (http://www.mathworks.com/help/stats/negative-binomial-distribution.html)
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = nbinpdf(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitlognpdf(tau,x,y)
% Fitting function for iterative fit to the sum of
% Lognormal Distributions
% (http://www.mathworks.com/help/stats/lognormal-distribution.html)
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = lognormal(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = lognormal(x,pos,wid)
% lognormal function.  pos=position; wid=half-width (both scalar)
% lognormal(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991  
g = exp(-(log(x/pos)/(0.01.*wid)) .^2);
% ----------------------------------------------------------------------
function err = fitsine(tau,x,y)
% Fitting function for iterative fit to the sum of
% sine waves (alpha test, NRFPT)
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = sine(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g=sine(x,f,phase) 
% Sine wave (alpha test)
g=sin(2*pi*f*(x+phase));
% ----------------------------------------------------------------------
function err = fitd1gauss(lambda,t,y)
%   Fitting functions for the first derivative of a Gaussian
%  T. C. O'Haver, 2014
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = d1gauss(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function y=d1gauss(x,p,w)
% First derivative of Gaussian (alpha test)
y=-(5.54518.*(x-p).*exp(-(2.77259.*(p-x).^2)./w^2))./w.^2;
y=y./max(y);
% ----------------------------------------------------------------------
function coeff = fitpolynomial(t,y,order)
coeff=polyfit(t,y,order);
% order=order
% coeff=coeff
% ----------------------------------------------------------------------
function y=polynomial(t,coeff)
y=polyval(coeff,t);
% ----------------------------------------------------------------------
function err = fitsegmented(lambda,t,y)
%   Fitting functions for articulated segmented linear
%  T. C. O'Haver, 2014
global LOGPLOT
breakpoints=[t(1) lambda max(t)];
z = segmented(t,y,breakpoints);
% lengthz=length(z);
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y);
end
% ----------------------------------------------------------------------
function yi=segmented(x,y,segs)
global PEAKHEIGHTS
clear yy
for n=1:length(segs)
  yind=val2ind(x,segs(n));
  yy(n)=y(yind(1));
end
yi=INTERP1(segs,yy,x);
PEAKHEIGHTS=segs;
% ----------------------------------------------------------------------
function err = fitlinslope(tau,x,y)
% Fitting function for iterative fit to linear function
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    z = (x.*tau(2*j-1)+tau(2*j))';
    A(:,j) = z./max(z);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function y=linslope(x,slope,intercept)
y=x.*slope+intercept;
% y=y./max(y);
% ----------------------------------------------------------------------
function b=iqr(a)
% b = IQR(a)  returns the interquartile range of the values in a.  For
%  vector input, b is the difference between the 75th and 25th percentiles
%  of a.  For matrix input, b is a row vector containing the interquartile
%  range of each column of a.
%  T. C. O'Haver, 2012
mina=min(a);
sizea=size(a);
NumCols=sizea(2);
for n=1:NumCols,b(:,n)=a(:,n)-mina(n);end
Sorteda=sort(b);
lx=length(Sorteda);
SecondQuartile=round(lx/4);
FourthQuartile=3*round(lx/4);
b=abs(Sorteda(FourthQuartile,:)-Sorteda(SecondQuartile,:));
% ----------------------------------------------------------------------
function err = fitmultiple(lambda,t,y,shapesvector,m)
% Fitting function for a multiple-shape band signal.
% The sequence of peak shapes are defined by the vector "shape".
% The vector "m" determines the shape of variable-shape peaks.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT coeff
numpeaks=round(length(lambda)/2);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    if shapesvector(j)==28,
        coeff=polyfit(t,y,m(j));
        A(:,j) = polyval(coeff,t);
    else
        A(:,j) = peakfunction(shapesvector(j),t,lambda(2*j-1),lambda(2*j),m(j))';
    end
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function p=peakfunction(shape,x,pos,wid,m,coeff)
% function that generates any of 20 peak types specified by number. 'shape'
% specifies the shape type of each peak in the signal: "peakshape" = 1-20.
% 1=Gaussian 2=Lorentzian, 3=logistic, 4=Pearson, 5=exponentionally
% broadened Gaussian; 9=exponential pulse, 10=up sigmoid,
% 13=Gaussian/Lorentzian blend; 14=BiGaussian, 15=Breit-Wigner-Fano (BWF) ,
% 18=exponentionally broadened Lorentzian; 19=alpha function; 20=Voigt
% profile; 21=triangular; 23=down sigmoid; 25=lognormal. "m" is required
% for variable-shape peaks only.
switch shape,
    case 1
        p=gaussian(x,pos,wid);
    case 2
        p=lorentzian(x,pos,wid);
    case 3
        p=logistic(x,pos,wid);
    case 4
        p=pearson(x,pos,wid,m);
    case 5
        p=expgaussian(x,pos,wid,m);
    case 6
        p=gaussian(x,pos,wid);
    case 7
        p=lorentzian(x,pos,wid);
    case 8
        p=expgaussian(x,pos,wid,m)';
    case 9
        p=exppulse(x,pos,wid);
    case 10
        p=upsigmoid(x,pos,wid);
    case 11
        p=gaussian(x,pos,wid);
    case 12
        p=lorentzian(x,pos,wid);
    case 13
        p=GL(x,pos,wid,m);
    case 14
        p=BiGaussian(x,pos,wid,m);
    case 15
        p=BWF(x,pos,wid,m);
    case 16
        p=gaussian(x,pos,wid);
    case 17
        p=lorentzian(x,pos,wid);
    case 18
        p=explorentzian(x,pos,wid,m)';
    case 19
        p=alphafunction(x,pos,wid);
    case 20
        p=voigt(x,pos,wid,m);
    case 21
        p=triangular(x,pos,wid);    
    case 23
        p=downsigmoid(x,pos,wid);
    case 25
        p=lognormal(x,pos,wid);
    case 26
        p=linslope(x,pos,wid);
    case 27
        p=d1gauss(x,pos,wid);
    case 28
        p=polynomial(x,coeff);
    otherwise
end % switch
