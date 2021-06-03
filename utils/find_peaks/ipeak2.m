function P=ipeak2(DataMatrix,PeakD,AmpT,SlopeT,SmoothW,FitW)
% ipeak2(DataMatrix)
%  Keyboard-operated Interactive Peak Finder for data in a single data
%  matrix "DataMatrix", with x values in row 1 and y values in row 2.
%  Experimental resolution enhancement function added June 2011.
%  Returns the peak table in P (Peak #, Position, Height, Width, Area)
%  http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
%  T. C. O'Haver (toh@umd.edu), Version 1.01, June 2011.  
%  
% EXAMPLE 1:   x=[0:.1:100];
%              y=(x.*sin(x)).^2;
%              ipeak2(x,y);  % Data in separate x and y vectors
% 
% EXAMPLE 2:   datamatrix=[x;y];  
%              ipeak2(datamatrix); % Data in two columns of a matrix
%
% EXAMPLE 3:   x=[0:.1:100];
%              y=5+5.*cos(x)+randn(size(x));
%              ipeak2(x,y,10);
%          or  ipeak2(datamatrix,10);
% The additional argument - 10 in this case - is an estimate of the peak
% density (maximum number of peaks that would fit into the data record).
% Small values detect fewer peaks; larger values > more peaks.
%
% EXAMPLE 4: ipeak2(datamatrix,0,.5,.0001,33,33);
%         or ipeak2(x,y,0,.5,.0001,33,33);
% As above, but specifies initial values of AmpT, SlopeT, SmoothW, FitW.
%
% Keyboard Controls:
% Pan signal left and right...Coarse pan: < and >   
%                             Fine pan: left and right cursor arrow keys
% Zoom in and out.............Coarse zoom: / and '   
%                             Fine zoom: up and down cursor arrow keys
% Adjust AmpThreshold.........A,Z  (A increases, Z decreases)
% Adjust SlopeThreshold.......S,X  (S increases, X decreases)
% Adjust SmoothWidth..........D,C  (D increases, C decreases)
% Adjust FitWidth.............F,V  (F increases, Z decreases)
% Baseline....................B, then click baseline at 8 points
% Restore original signal.....G
% Print peak table............P  Prints Peak #, Position, Height, Width
% Print keyboard commands.....K  Prints this list
% Print findpeaks arguments...Q   (AmpT, SlopeT, SmoothW, FitW)
% Print report................R  Prints Peak table and parameters
% Peak labels ON/OFF......... L  Displays peak parameters in upper window.
%
% For large data sets, to view only a portion of the data over a 
% specified x-axis range, you can type 
% n1=val2ind(x,x1);n2=val2ind(x,x2);ipeak([x(n1:n2)' y(n1:n2)'])
% where x1 and x2 are the end limits of the x-axis values displayed.
global X
global Y
global xo
global dx
global SlopeThreshold 
global AmpThreshold  
global SmoothWidth
global FitWidth
global PeakLabels
global dsmooth
global REfactor1
global REfactor2
global P
% close
format short g
format compact
warning off all

switch nargin
    % 'nargin' is the number of arguments
    case 1
      % Assumne that the argument must be a matrix of data.
      % If DataMatrix is in the wrong transposition, fix it.
      datasize=size(DataMatrix);
      if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
      % Assign arguments to internal global variables
      X=DataMatrix(:,1); % Split matrix argument 
      Y=DataMatrix(:,2);
      % Calculate default values of peak detection parameters
      PeakDensity=20;   
      % Estimate approx number of points in a peak half-width
      WidthPoints=length(Y)/PeakDensity;  
      SlopeThreshold=WidthPoints^-2;  
      AmpThreshold=abs(min(Y)+0.05*(max(Y)-min(Y))); 
      SmoothWidth=round(WidthPoints/3);  
      FitWidth=round(WidthPoints/3); 
  case 2
      % Two arguments, might be separate x and y data vectors, 
      % or one data matrix and a peak density estimate.
      if isscalar(PeakD)
         % Must be one data matrix and a peak density estimate.
         % If DataMatrix is in the wrong transposition, fix it.
         datasize=size(DataMatrix);
         if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
         % Calculate values of peak detection parameters
         % arguments based on the peak density, PeakD
         size(DataMatrix)
         PeakDensity=PeakD;    
         % Estimate approx number of points in a peak half-width
         WidthPoints=length(Y)/PeakDensity;  
         SlopeThreshold=WidthPoints^-2;  
         AmpThreshold=abs(min(Y)+0.02*(max(Y)-min(Y)));  
         SmoothWidth=round(WidthPoints/3);  
         FitWidth=round(WidthPoints/3); 
      else % if not isscalar
        % Must be separate x and y data vectors
        X=DataMatrix;
        Y=PeakD;
        PeakDensity=20;   
        % Estimate approx number of points in a peak half-width
        WidthPoints=length(Y)/PeakDensity;  
        SlopeThreshold=WidthPoints^-2;  
        AmpThreshold=abs(min(Y)+0.05*(max(Y)-min(Y))); 
        SmoothWidth=round(WidthPoints/3);  
        FitWidth=round(WidthPoints/3); 
      end  % if isscalar
    case 3
      % Must be  be separate x and y data vectors plus  a peak density estimate. 
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
    case 6 
      % Must be one data matrix and all parameters specified in arguments 
      % If DataMatrix is in the wrong transposition, fix it.
      datasize=size(DataMatrix);
      if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
      X=DataMatrix(:,1); % Split matrix argument 
      Y=DataMatrix(:,2);
      SlopeThreshold=SlopeT;
      AmpThreshold=AmpT;
      SmoothWidth=SmoothW;
      FitWidth=FitW; 
    case 7
      % Must be separate x and y data vectors and all parameters specified in arguments 
      X=DataMatrix;
      Y=PeakD;
      SlopeThreshold=SlopeT;
      AmpThreshold=AmpT;
      SmoothWidth=SmoothW;
      FitWidth=FitW;  
  otherwise
      disp('Invalid number of arguments')
      disp('Expected forms are:')
      disp('ipeak(x,y);  % Data in separate x and y vectors')
      disp('ipeak(datamatrix); % Data in two columns of datamatrix')
      disp('ipeak(x,y,10); or ipeak(datamatrix,10);  specifying peak density')
      disp('ipeak(x,y,0,.5,.0001,33,33);  specifying peak density, AmpT, SlopeT, SmoothW, FitW')
      disp('ipeak(datamatrix,0,.5,.0001,33,33);')
      return
end % switch nargin

if FitWidth<3,FitWidth=3;end   % Keep FitWidth above 2 
xo=length(Y)/2; % Initial Pan setting
dx=length(Y)/4; % Initial Zoom setting
PeakLabels=0; % Peak numbers only, no parameter labels, in upper window
dsmooth=SmoothWidth; % initial dsmooth setting for resolution enhancement
REfactor1=dsmooth^2/25; % factor1 for resolution enhancement (Gaussian)
REfactor2=dsmooth^4/833; % factor2 for resolution enhancement  (Gaussian)
% REfactor1=dsmooth^2/6; % factor1 for resolution enhancement (Lorentzian)
% REfactor2=dsmooth^4/700; % factor2 for resolution enhancement  (Lorentzian) 

% Plot the signal 
RedrawSignal(X,Y,xo,dx);

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
global X
global Y
global xx
global yy
global xo
global dx
global SlopeThreshold 
global AmpThreshold  
global SmoothWidth
global FitWidth
global PeakLabels
global SavedSignal
global oldAmpThreshold
global dsmooth
global REfactor1
global REfactor2
global P
key=get(gcf,'CurrentCharacter');
if ischar(key),
  switch double(key),
    case 29
        % Pans one point down when left arrow pressed.
          xo=xo-1;
          if xo<1,xo=1;end
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 28
        % Pans one point up when right arrow pressed.
        ly=length(Y);
        xo=xo+1;
        if xo>ly,xo=ly;end
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 46
        % Pans 2% down when < key pressed.
        ly=length(Y);
        xo=xo-ly/50;
        if xo<1,xo=1;end
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 44
        % Pans 2% up when > key pressed.
        ly=length(Y);
        xo=xo+ly/50;
        if xo>ly,xo=ly;end
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 31
        % Zooms one point up when up arrow pressed.
        dx=dx+2;
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 30
        % Zooms one point down when down arrow pressed.
        dx=dx-2;
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 47
        % Zooms 2% up when / pressed.
        ly=length(Y);
        dx=dx+ly/50;
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 39
        % Zooms 2% down when ' pressed.
        ly=length(Y);
        dx=dx-ly/50;
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 98
        % When 'b' key is pressed, user clicks graph 
        % to enter background points, then graph re-drawn.
        SavedSignal=Y;
        oldAmpThreshold=AmpThreshold;
        BaselinePoints=8;  % Change as you wish <============
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
        % Estimate initial value of AmpThreshold
        AmpThreshold=min(Y)+0.02*(max(Y)-min(Y));  
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
      case 103
          % When 'g' key is pressed, restores signal background and AmpThreshold. 
          Y=SavedSignal;
          AmpThreshold=oldAmpThreshold;
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
      case 97
        % When 'a' key is pressed, increases "AmpThreshold" by 10%
        AmpThreshold=abs(AmpThreshold+.1*AmpThreshold);
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 122
        % When 'z' key is pressed, decreases "AmpThreshold" by 10%
        AmpThreshold=AmpThreshold-.1*AmpThreshold;
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 115 % When 's' key is pressed, increases "SlopeThreshold" by 10%
         SlopeThreshold=SlopeThreshold+.1*SlopeThreshold;
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 120 % When 'x' key is pressed, decreases "SlopeThreshold" by 10%
         SlopeThreshold=SlopeThreshold-.1*SlopeThreshold;
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 100
        % When 'd' key is pressed, increases "SmoothWidth" by 1
        SmoothWidth=SmoothWidth+1;
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 99
        % When 'c' key is pressed, decreases "SmoothWidth" by 1
        SmoothWidth=SmoothWidth-1;
        if SmoothWidth<1, SmoothWidth=1;end
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 102
        % When 'f' key is pressed, increases "FitWidth" by 1
        FitWidth=FitWidth+1;
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 118
        % When 'v' key is pressed, decreases "FitWidth" by 1
         FitWidth=FitWidth-1;
         if FitWidth<0, FitWidth=0;end
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 104
        % When 'h' key is pressed, increases REfactors  1 and 2
        REfactor1=REfactor1+.1*REfactor1;
        REfactor2=REfactor2+.1*REfactor2;
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 110
        % When 'n' key is pressed, decreases REfactors 1 and 2
         REfactor1=REfactor1-.1*REfactor1;
         if REfactor1<0, REfactor1=0;end
         REfactor2=REfactor2-.1*REfactor2;
         if REfactor2<0, REfactor2=0;end
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);       
    case 106
        % When 'j' key is pressed, increases "dsmooth"
        dsmooth=dsmooth+2;
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
    case 109
        % When 'm' key is pressed, decreases "dsmooth"
         dsmooth=dsmooth-2;
         if dsmooth<0, dsmooth=0;end
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);       
    case 121
        % When 'y' key is pressed, cancels resolution enhancement
        REfactor1=1/25;
        REfactor2=1/800;
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);   
      case 114
        % When 'r' key is pressed, prints
        disp('--------------------------------------------------------')
        disp(['Amplitude Threshold (AmpT) = ' num2str(AmpThreshold) ] )
        disp(['Slope Threshold (SlopeT) = ' num2str(SlopeThreshold) ] )
        disp(['Smooth Width (SmoothW) = ' num2str(SmoothWidth) ] )
        disp(['Resolution Enhancement factor 1(REf) = ' num2str(REfactor1) ] )
        disp(['RE factor2 = ' num2str(REfactor2) ] )
        disp(['RE smooth (REs) = ' num2str(dsmooth) ] )
        disp('         Peak #    Position      Height       Width        Area')
        disp(P)
    case 112
        % When 'p' key is pressed, prints out peak table 
        disp('         Peak #    Position      Height       Width        Area')
        disp(P)
    case 111
        % When 'o' key is pressed, the signal Y is saved to disk as
        % SignalOutput.mat
        output=enhance(Y,REfactor1,REfactor2,dsmooth);
        save SignalOutput output
    case 107
        % When 'k' key is pressed, prints out table of keyboard commands
        disp('KEYBOARD CONTROLS:')
        disp(' Pan signal left and right...Coarse pan: < and >')   
        disp('                             Fine pan: left and right cursor arrows')
        disp(' Zoom in and out.............Coarse zoom: / and "  ') 
        disp('                             Fine zoom: up and down cursor arrows')
        disp(' Adjust AmpThreshold.........A,Z (Larger values ignore short peaks)')
        disp(' Adjust SlopeThreshold.......S,X (Larger values ignore broad peaks)')
        disp(' Adjust SmoothWidth..........D,C (Larger values ignore sharp peaks}')
        disp(' Adjust FitWidth.............F,V (Adjust to cover just top part of peaks')
        disp(' Adjust REf..................H,N  Resolution enhancement strength')
        disp(' Adjust REs..................J,M  Resolution enhancement smoothing')
        disp(' Cancel resolution enhancement...Y  ')
        disp(' Baseline subtraction........B, then click baseline at 8 points')
        disp(' Restore original signal.....G  to cancel previous background subtraction')
        disp(' Print report................R  Prints Peak table and parameters')
        disp(' Print peak table............P  Peak #, Position, Height, Width, Area')
        disp(' Print keyboard commands.....K  Prints this list')
        disp(' Print findpeaks arguments...Q  (AmpT, SlopeT, SmoothW, FitW)')
        disp(' Peak labels ON/OFF..........L  Displays peak parameters on upper graph')
    case 113
        % When 'Q' is pressed, prints peak detection parameters on a single line     
        disp([ num2str(AmpThreshold) ',' num2str(SlopeThreshold)  ',' num2str(SmoothWidth)  ',' num2str(FitWidth) ])
    case 108
        % When 'L' is pressed, toggles on/off peak labels in upper panel
        if PeakLabels==0,
            PeakLabels=1;
        else
            PeakLabels=0;
        end
        [xx,yy]=RedrawSignal(X,enhance(Y,REfactor1,REfactor2,dsmooth),xo,dx);       
    otherwise  
       UnassignedKey=double(key)  % Display key code of unassigned keys
       disp('Press k to print out list of keyboard commands')
   end % switch
end % if
% ----------------------------------------------------------------------    
function [xx,yy]=RedrawSignal(X,Y,xo,dx)
% Plots the entire signal (X,Y) in the lower half of the plot window and an
% isolated segment (xx,yy) in the upper half, controlled by Pan and Zoom
% keys.
global SlopeThreshold 
global AmpThreshold  
global SmoothWidth
global FitWidth
global PeakLabels
global dsmooth
global REfactor1
global REfactor2
global P
Startx=round(xo-(dx/2));
Endx=round(xo+(dx/2)-1);
if Endx>length(Y),
    Endx=length(Y);
end
if Startx<1,
     Startx=1;
end
PlotRange=[Startx:Endx];
if PlotRange<5, PlotRange=[xo:xo+5];end
xx=X(PlotRange);
yy=Y(PlotRange); 
hold off
% clf
% Plots isolated segment (xx,yy) in the upper half
figure(1);subplot(2,1,1);plot(xx,yy,'g.'); 
hold on
title('Press < > to pan and ? " to zoom. Press K for keyboard commands')
axis([X(Startx) X(Endx) min(yy) max(yy)]);

% Bottom half of the figure shows full signal
figure(1);subplot(2,1,2);cla
plot(X,Y,'b')  % Graph the signal
axis([X(1) X(length(X)) min(Y) max(Y)]); % Update plot
P=findpeaks(X,Y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth);
title(' AmpT: A/Z  SlopeT: S/X  SmoothWidth: D/C  FitWidth: F/V  REf: H/N   REs:J/M')
xlabel(['    AmpT: ' num2str(AmpThreshold) '    SlopeT: ' num2str(SlopeThreshold) '   SmoothW: ' num2str(SmoothWidth) '   FitW: ' num2str(FitWidth) '   REf: ' num2str(REfactor1) '    REs: ' num2str(dsmooth) ])
text(P(:,2),P(:,3),num2str(P(:,1)),'FontSize',12)  % Number the peaks found on the lower grap
hold on
% Mark the zoom range on the full signal with two magenta dotted vertical lines
plot([min(xx) min(xx)],[min(Y) max(Y)],'m--')
plot([max(xx) max(xx)],[min(Y) max(Y)],'m--') 
% Number the peaks found on the upper graph
subplot(2,1,1);
if PeakLabels==1,
   % Label the peaks on the upper graph with number, position, height, and
   % width
   topaxis=axis;
   yrange=topaxis(4)-topaxis(3);
   pos1=.1*yrange;
   pos2=.2*yrange;
   pos3=.3*yrange;
   pos4=.4*yrange;
   pos5=.5*yrange;
   PP=findpeaks(xx,yy,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth);
   text(P(:,2),P(:,3)-pos1,num2str(P(:,1)),'FontSize',12)
   text(PP(:,2),PP(:,3)-pos2,num2str(PP(:,2)))
   text(PP(:,2),PP(:,3)-pos3,num2str(PP(:,3)))
   text(PP(:,2),PP(:,3)-pos4,num2str(PP(:,4)))
   text(PP(:,2),PP(:,3)-pos5,num2str(PP(:,5)))
else   
   topaxis=axis;
   yrange=topaxis(4)-topaxis(3);
   pos1=.1*yrange; 
   % Number the peaks on the upper graph
   text(P(:,2),P(:,3)-pos1,num2str(P(:,1)),'FontSize',12) 
end
markpeaks(xx,yy,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth);
hold off
% ----------------------------------------------------------------------
function markpeaks(xx,yy,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth)
% Draws best-fit line through data points near peak and computes
% peak position, amplitude, and width
% PP=[0 0 0 0];
PP=findpeaks(xx,yy,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth);
hold on
if PP(1)>0,
   sizePP=size(PP);
   lengthPP=sizePP(1);
  for PeakNumber=1:lengthPP,
    FitWidth=round(real(FitWidth));
    subplot(2,1,1);
    if PeakNumber>lengthPP,PeakNumber=lengthPP;end
    n1=val2ind(xx,PP(PeakNumber,2))-round(FitWidth/2);
    n2=val2ind(xx,PP(PeakNumber,2))+round(FitWidth/2);
    if n1<1, n1=1;end
    if n2>length(yy), n2=length(yy);end
    PlotRange=[n1:n2];
    xxx=xx(PlotRange);
    yyy=yy(PlotRange);
    % Fit parabola to log10 of sub-group
    logyyy=log(abs(yyy));
    if size(xxx)~=size(logyyy),logyyy=logyyy';end
    [coef,S,MU]=polyfit(xxx,logyyy,2);  
    c1=coef(3);c2=coef(2);c3=coef(1);
    % Compute peak position and height or fitted parabola
    PeakX=-((MU(2).*c2/(2*c3))-MU(1));  
    PeakY=exp(c1-c3*(c2/(2*c3))^2);
    MeasuredWidth=norm(MU(2).*2.35703/(sqrt(2)*sqrt(-1*c3)));
    subplot(2,1,1);
    plot(xxx,PeakY.*gaussian(xxx,PeakX,MeasuredWidth),'k');
  end  
end
% ----------------------------------------------------------------------
function P=findpeaks(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup)
% function P=findpeaks(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup)
% Function to locate the positive peaks in a noisy x-y data
% set.  Detects peaks by looking for downward zero-crossings
% in the first derivative that exceed SlopeThreshold.
% Returns list (P) containing peak number and the position, 
% height, width, and area of each peak. SlopeThreshold,
% AmpThreshold, and smoothwidth control sensitivity
% Higher values will neglect smaller features. Peakgroup
% is the number of points around the "top part" of the peak.
% T. C. O'Haver, 1995.  Version 3; revised June 7, 2011
smoothwidth=round(smoothwidth);
peakgroup=round(peakgroup);
d=fastsmooth(deriv(y),smoothwidth);
n=round(peakgroup/2+1);
P=[0 0 0 0 0];
vectorlength=length(y);
peak=1;
AmpTest=AmpThreshold;
for j=smoothwidth:length(y)-smoothwidth,
   if sign(d(j)) > sign (d(j+1)), % Detects zero-crossing
       % if slope of derivative is larger than SlopeThreshold
     if d(j)-d(j+1) > SlopeThreshold*y(j), 
        if y(j) > AmpTest,  % if height of peak is larger than AmpThreshold
          for k=1:peakgroup, % Create sub-group of points near peak
              groupindex=j+k-n+1;
              if groupindex<1, groupindex=1;end
              if groupindex>vectorlength, groupindex=vectorlength;end
            xx(k)=x(groupindex);yy(k)=y(groupindex);
          end% Fit parabola to log10 of sub-group with centering and scaling
          [coef,S,MU]=polyfit(xx,log(abs(yy)),2);  
          c1=coef(3);c2=coef(2);c3=coef(1);
          % Compute peak position and height of fitted parabola
          PeakX=-((MU(2).*c2/(2*c3))-MU(1));   
          PeakY=exp(c1-c3*(c2/(2*c3))^2);
          MeasuredWidth=norm(MU(2).*2.35703/(sqrt(2)*sqrt(-1*c3)));
          % if the peak is too narrow for least-squares technique to work
          % well, just use the max value of y in the sub-group of points near peak.
          if peakgroup<7,
             PeakY=max(yy);
             pindex=val2ind(yy,PeakY);
             PeakX=xx(pindex(1));
          end         
         % Construct matrix P. One row for each peak detected, 
          % containing the peak number, peak position (x-value),  
          % height (y-value), width, and area (height x width).
          P(peak,:) = [round(peak) PeakX PeakY MeasuredWidth PeakY*MeasuredWidth];
          peak=peak+1;
        end
      end
   end
end
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
function g = gaussian(x,pos,wid)
% gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
% X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6006.*wid)) .^2);
% ----------------------------------------------------------------------
function Enhancedsignal=enhance(signal,factor1,factor2,SmoothWidth)
% Resolution enhancement function by even derivative method. The
% arguments factor1 and factor 2 are 2nd and 4th derivative weighting
% factors. Larger values of factor1 and factor2 will reduce the 
% peak width but will cause artifacts in the baseline near 
% the peak.  Adjust the factors for the the best compromise. 
% Use minimum smooth width needed to reduce excess noise. 
d2=secderiv(signal');  % Computes second derivative
d4=secderiv(d2);   % Computes fourth derivative
Enhancedsignal = signal'-factor1.*fastsmooth(d2,SmoothWidth,2)+...
factor2.*fastsmooth(fastsmooth(fastsmooth(d4,SmoothWidth,2),SmoothWidth,2),SmoothWidth,2);
% ----------------------------------------------------------------------
function s=tsmooth(Y,w)
%  tsmooth(Y,w) smooths vector Y by a triangular function of halfwidth w
%  T. C. O'Haver, 1988.
v=conv(boxcar(w),boxcar(w));
S=conv(Y,v);
startpoint=(length(v) + 1)/2;
endpoint=length(Y)+startpoint-1;
s=S(startpoint:endpoint) ./ sum(v);
% ----------------------------------------------------------------------
function s=boxcar(w)
%  boxcar(w) = Rectangular function of width w
%  T. C. O'Haver, 1988.
s=ones(1,w);
% ----------------------------------------------------------------------
function d=secderiv(a)
% Second derivative of vector using 3-point central difference.
%  T. C. O'Haver, 2006.
n=length(a);
for j = 2:n-1;
  d(j)=a(j+1) - 2.*a(j) + a(j-1);
end
d(1)=d(2);
d(n)=d(n-1);
% ----------------------------------------------------------------------
function d=deriv(a)
% First derivative of vector using 2-point central difference.
%  T. C. O'Haver, 1988.
n=length(a);
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1;
  d(j)=(a(j+1)-a(j-1)) ./ 2;
end
% ----------------------------------------------------------------------
function SmoothY=fastsmooth(Y,w,type,ends)
% fastbsmooth(Y,w,type,ends) smooths vector Y with smooth 
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
