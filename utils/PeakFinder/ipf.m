function ipf(arg1,arg2,arg3,arg4)
% Keyboard-operated Interactive Peak Fitter for data in separate x,y
% vectors, or in a single y vector, or in a 2xn or nx2 data matrix with
% x values in row/column 1 and y values in row/column 2 (e.g. [x y]).
% Press K key for list of keyboard commands.
% Version 12.3: October, 2017, Changed IQR in the bootstrap method to
% IQR/1.34896, which estimates the RSD without outliers for a normal
% distribution.
% Syntax variations:
%  One input argument:
%   ipf(y) or ipf([x;y]) or ipf([x;y]');
%  Two input arguments:
%   ipf(x,y) or ipf([x;y],center) or ipf([x;y]',center);
%  Three input arguments:
%   ipf(x,y,center) or ipf(y,center,window) or 
%   ipf([x;y],center,window) or ipf([x;y]',center,window);
%  Four input arguments:
%   ipf(x,y,center,window);
% See http://terpconnect.umd.edu/~toh/spectrum/ifpinstructions.html
% See http://terpconnect.umd.edu/~toh/spectrum/CurveFittingC.html
% See http://www.wam.umd.edu/~toh/spectrum/InteractivePeakFitter.htm
% T. C. O'Haver (toh@umd.edu).
%
% Example 1: 
% x=[0:.005:1];y=humps(x);ipf(x,y)
% Example 2: 
% x=[-10:.1:10];y=exp(-(x).^2);ipf(x,y)
% Example 3:
% x=[0:.1:10];y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.1*randn(size(x));ipf(x,y)
% Example 4:
% x=[0:.01:18];
% y=exp(-(x-4).^2)+exp(-(x-9).^2)+exp(-(x-12).^2)+exp(-(x-13.7).^2);
% ipf(x,y,11,10);
% Example 5
% x=[-.1:.005:1.2];y=humps(x);ipf(x,y,0.6,1.2)
% Example 6
% x=[-5:.02:5];y=100.*x+100.*exp(-(x).^2)+randn(size(x));ipf(x,y)
%
% Keyboard Controls:
% Pan signal left and right...Coarse pan: < and >
%                             Fine pan: left and right cursor arrow keys
%                             Nudge: [ ]
% Zoom in and out.............Coarse zoom: / and '
%                             Fine zoom: up and down cursor arrow keys
% Resets pan and zoom.........ESC
% Select entire signal........Crtl-A
% Select # of peaks...........Number keys 1-9, or press 0 key to enter manually
% Complete peak shape menu... - (minus or hyphen), then type number and Enter
% Select peak shape...........g Gaussian
%                             G (Shift-G) Fixed-width Gaussian
%                             Shift-P  Fixed-position Gaussians
%                             h Equal-width Gaussian
%                             H (Shift-H) bifurcated Gaussian (a,z keys adjust shape)
%                             l Lorentzian
%                             ; Equal-width Lorentzian
%                             L (Shift-L) Fixed-width Lorentzian
%                             Shift-[  Fixed-position Lorentzians
%                             Shift-B  Breit-Wigner-Fano (a,z keys adjust Fano factor)
%                             o Lognormal distribution (See "S" for logistic function)
%                             p Pearson (use a,z keys to adjust shape)
%                             e exponentially-broadened Gaussian
%                             E (Shift-E) Exponential-broadened Lorentzians
%                             j exponential-broadened equal-width Gaussian
%                                 (a,z keys adjust broadening)
%                             u exponential pUlse: y=exp(-tau1.*x).*(1-exp(-tau2.*x))
%                             U (Shift-U) Alpha: y=(x-tau2)./tau1.*exp(1-(x-tau2)./tau1);
%                             s Up Sigmoid (logistic function): y=.5+.5*erf((x-tau1)/sqrt(2*tau2))                                                       
%                             Shift-D Down Sigmoid (logistic function):  y=.5-.5*erf((x-tau1)/sqrt(2*tau2))                                                      
%                             ~ Gauss/Lorentz blend (a/z keys adjust shape)
%                             V (Shift-V) Voigt profile (a/z adjusts shape)
%                             Shift-T  Triangular
% Fit.........................f  perform single Fit from another start point.
% Select autozero mode........t  selects none, linear, quadratic, or flat baseline mode
% Monopolar/bipolar mode......= + Flips between + peaks only and +/- peak mode
% Toggle log y mode OFF/ON....m  Plot linear or log Y axis on lower graph
% 2-point Baseline............b, then click left and right baseline
% Set manual baseline.........Backspace, then click baseline at multiple points
% Restore original baseline...\  to cancel previous background subtraction
% Click start positions.......c, click on peak position for each peak
% Type in start vector........C (Shift-C) Type or Paste start vector [p1 w1 p2 w2 ...]
% Enter value of 'extra'......Shift-x, type value, press Enter.
% Adjust 'extra' up/down......a,z: 5% change; upper case A,Z: 0.5% change.
% Print parameters & results..q
% Print fit results only......r
% Plot signal in figure 2.....y
% Print out model data table..d
% Refine fit..................x  Takes best of 10 trial fits (change in line 224)
% Print peakfit function......w  with all input arguments
% Compute bootstrap stats.....v  Estimates standard deViations of parameters.
% Test effect of Noise........n  Test effect of Noise by fitting subset of points
% Save Figure as png file.....Shift-S  Saves as Figure1.png, Figure2.png....
% Display current settings....Shift-? displays table of current values
% Prints list of commands.....k
% Switch to iPeak.............Shift-Ctrl-P transfer current signal to iPeak.m
% Switch to iSignal...........Shift-Ctrl-S transfer current signal to iSignal.m
% Fit polynomial to segment...Shift-o  Asks for polynomial order
% Enter minimun width.........Shift-W  
% Enter saturation maximum....Shift-M
%
% Copyright (c) 2015, Thomas C. O'Haver
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
% 
global X Y xx yy xo dx NumPeaks NumTrials Shape shapesvector AA PEAKHEIGHTS xxx FigNum
global start extra delta AUTOZERO SavedSignal logplot FIXEDPARAMETERS BIPOLAR
global CLIPHEIGHT MINWIDTH
%
format short g
format compact
warning off all
warning OFF BACKTRACE
warning OFF VERBOSE

% process input arguments
switch nargin % Process arguments
    % 'nargin' is the number of arguments
    case 1  % One argument only
        % Might be ipf(DataMatrix) ot ipf(Y-vector)
        % If data is in the wrong transposition, fix it.
        datasize=size(arg1);
        if datasize(1)<datasize(2),arg1=arg1';end
        datasize=size(arg1);
        if datasize(2)==1, %  Must be ipf(Y-vector)
            X=1:length(arg1); % Create an independent variable vector
            Y=arg1;
        else
            % Must be ipf(DataMatrix)
            X=arg1(:,1); % Split matrix argument
            Y=arg1(:,2);
        end
        xo=length(Y)/2; % Initial Pan setting
        dx=length(Y)/4; % Initial Zoom setting
    case 2
        % Two arguments, might be separate x and y data vectors,
        % or one adata matrix and a peak density estimate.
        if isscalar(arg2) % if second argument is scalar
            % Must be ipf(DataMatrix,center)
            % If DataMatrix is in the wrong transposition, fix it.
            datasize=size(arg1);
            if datasize(1)<datasize(2),arg1=arg1';end
            X=arg1(:,1); % Split matrix argument
            Y=arg1(:,2);
            xo=val2ind(X,arg2);
        else % if second argument is not scalar
            % Must be ipf(x,y)
            xdatasize=size(arg1);
            if xdatasize(1)<xdatasize(2),arg1=arg1';end
            X=arg1;  % First argument is X
            xdatasize=size(arg2);
            if xdatasize(1)<xdatasize(2),arg2=arg2';end
            Y=arg2; % Second argument is Y
            xo=length(Y)/2; %  % Default initial zoom setting
        end  % if isscalar
        dx=length(Y)/4; %  % Default initial zoom setting
    case 3
        % Might be ipf(DataMatrix,center,window) or ipf(x,y,center)
        if isscalar(arg2) % if second argument is scalar
            % Must be ipf(DataMatrix,center,window)
            % If DataMatrix is in the wrong transposition, fix it.
            datasize=size(arg1);
            if datasize(1)<datasize(2),arg1=arg1';'Transpose matrix';end
            datasize=size(arg1);
            if datasize(2)==2, % ipf(DataMatrix,center,window)
                X=arg1(:,1); % Split matrix argument
                Y=arg1(:,2);
                xo=val2ind(X,arg2); % Second argument is center
                dx=val2ind(X,arg3); % Third argument is window
            else % ipf(y,center,window)
                X=1:length(arg1); % Create an independent variable vector
                Y=arg1;
                xo=val2ind(X,arg2); % Second argument is center
                dx=val2ind(X,arg3); % Third argument is window
            end
        else % if second argument is not isscalar
            % Must be ipf(x,y,center)
            xdatasize=size(arg1);
            if xdatasize(1)<xdatasize(2),arg1=arg1';end
            X=arg1;  % First argument is X
            xdatasize=size(arg2);
            if xdatasize(1)<xdatasize(2),arg2=arg2';end
            Y=arg2; % Second argument is Y
            xo=val2ind(X,arg3); % third argument is center
            dx=length(Y)/4; % Default initial zoom setting
        end  % if isscalar
    case 4   % Must be ipf(x,y,center,window)
        % 'case 4   % Must be ipf(x,y,center,window) = ipf(DataMatrix,center,window,NumPeaks'
        xdatasize=size(arg1);
        if xdatasize(1)<xdatasize(2),arg1=arg1';end
        X=arg1;  % First argument is X
        xdatasize=size(X);
        if xdatasize(1)<xdatasize(2),X=X';end
        Y=arg2; % Second argument is Y
        xo=val2ind(X,arg3); % Third argument is center
        dx=val2ind(X,arg4); % Fourth argument is window
end

% Adjust X and Y vector shape to 1 x n (rather than n x 1)
X=reshape(X,1,length(X));
Y=reshape(Y,1,length(Y));
% If necessary, flip the data vectors so that X increases
if X(1)>X(length(X)),
    disp('X-axis flipped.')
    X=fliplr(X);
    Y=fliplr(Y);
end

% Set initial values of parameters
NumPeaks=1; % Initial Number of peaks in model
NumTrials=10; % Number of repeat fits when X key is pressed (CHANGE AS DESIRED)
Shape=1; % Initial Shape of the peaks (1=Gaussian, 2=Lorentzian, etc)
shapesvector=Shape; 
extra=length(X)/100.*ones(size(Shape));
if Shape==20,extra=1;end
FigNum=1;
delta=0; % Initial Random change in initial start guesses for each re-fit
PEAKHEIGHTS=zeros(NumPeaks,1);
xxx=zeros(1,200);
AA=zeros(NumPeaks,200);
AUTOZERO=0; % Default to no baseline correction Press T to select mode.
logplot=0; % Default to linear mode. Press m to change to log.
SavedSignal=Y;
FIXEDPARAMETERS=0; % Used only for fixed position or width shapes.
BIPOLAR=0; % Default to monopolar (positive peaks only) mode.
CLIPHEIGHT=max(Y);
MINWIDTH=X(2)-X(1);

% Plot the signal and its fitted components and residuals
[xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
% if start==[];start=length(Y)/2;end
FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
[xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);

% Attaches KeyPress test function to the figure.
set(gcf,'KeyPressFcn',@ReadKey)
uicontrol('Style','text')
% end of outer function
% ----------------------------SUBFUNCTIONS--------------------------------
function ReadKey(obj,eventdata)
% Interprets key presses from the Figure window.
global X Y xx yy xo dx start FitResults NumPeaks NumTrials Shape residual
global delta AA xxx PEAKHEIGHTS AUTOZERO extra MeanFitError SavedSignal
global  FIXEDPARAMETERS logplot et FigNum BIPOLAR shapesvector MINWIDTH CLIPHEIGHT
% (et=elapsed time)
% When a key is pressed, interprets key and calls corresponding function.
% Note: If you don't like my key assignments, you can change the numbers
% in the case statements here to re-assign that function to any other key.
NumTrialsBoot=1;
key=get(gcf,'CurrentCharacter');
if isscalar(key),
    ly=length(Y);
    switch double(key),
        case 28
            xo=xo-dx/20;
            if xo<1,xo=1;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 29
            xo=xo+dx/20;
            if xo>ly,xo=ly;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 91
            % Nudge down 1 point when [ pressed.
            xo=xo-1;
            if xo<1,xo=1;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 93
            % Nudge up 1 point when ] pressed.
            xo=xo+1;
            if xo>ly,xo=ly;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 44
            % Pans down when > key pressed.
            xo=xo-dx;
            if xo<1,xo=1;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 46
            % Pans up when < key pressed.
            xo=xo+dx;
            if xo>ly,xo=ly;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 31
            dx=dx+dx/10;
            if dx>ly,dx=ly;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 30
            dx=dx-dx/10;
            if dx<2,dx=2;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 47
            % Zooms 2% out when / pressed.
            dx=dx+ly/50;
            if dx>ly,dx=ly;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 39
            % Zooms 2% in when ' pressed.
            dx=dx-ly/50;
            if dx<2,dx=2;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 1 % Ctrl-A selects entire signal
            xo=length(Y)/2; % Initial Pan setting           
            dx=length(Y); % Initial Zoom setting
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);            
        case 109
            % When 'M' is pressed, toggles on/off log plot mode
            if logplot==0,
                logplot=1;
                SavedSignal=Y;
                Y=abs(Y);
            else
                logplot=0;
                Y=SavedSignal;
            end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 27 % When 'ESC' key is pressed, resets pan and zoom
            xo=length(Y)/2; % Initial Pan setting
            dx=length(Y)/4; % Initial Zoom setting
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 102
            % When 'f' key is pressed, tweaks start values, computes and
            % plots fit.
            center=(max(xx)+min(xx))/2;
            window=max(xx)-min(xx);
            startnow=start;
            delta=(max(xx)-min(xx))/100;
            for k=1:2*NumPeaks,
                startnow(k)=start(k)+(rand-.5)*delta;
            end
            NumTrials=1;
            [FitResults,MeanFitError]=peakfit([xx;yy],center,window,NumPeaks,shapesvector,extra,NumTrials,startnow,AUTOZERO,FIXEDPARAMETERS,1,BIPOLAR,MINWIDTH,delta);
            % [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,del
            % ta,startnow,extra);
           start=startnow;
        case 99
            % When 'c' key is pressed, user clicks graph to enter start positons,
            % then fit is computed and graph rte-drawn.
            % Acquire first-guess peak positions from user mouse pointer
            figure(1);subplot(2,1,1);xlabel('Click on the estimated positions of each proposed component peak.')
            [clickX,clickY] = ginput(NumPeaks);
            % Create a "start" vector using these peak positions, with peak
            % widths
            n=max(xx)-min(xx);
            width=n/(5*NumPeaks);
            start=[];
            for k=1:NumPeaks,
                start=[start clickX(k) width];
            end
            if Shape==11||12,
                fixedstart=[];
                for pk=1:NumPeaks,
                    fixedstart(pk)=start(2*pk-1);
                end
                peakfit([xx;yy],0,0,NumPeaks,shapesvector,extra,1,start,AUTOZERO,FIXEDPARAMETERS,1,BIPOLAR,MINWIDTH,delta);
            end
            [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
        case 67 % Shift-C allows keyboard entry of start vector
             disp('The current start vector is displayed below, in the order')
             disp('[position1 width1 position2 width2 ...]; you can Copy and Paste')
             disp('into the input prompt, edit desired values, and press Enter.')
             disp(['[' num2str(start) ']' ])
             startinput=input('Type [start vector] or press Enter to keep unchanged: ');
            if isempty(startinput),
            else
                start=startinput;
            end
             [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
        case 81 % Shift-Q
            disp(['First guess vector=' num2str(start)])
        case 98
            % When 'b' key is pressed, user clicks graph before and after peaks
            % to enter background points, then fit is computed and graph re-drawn.
            figure(1);subplot(2,1,1);xlabel('Click on the baseline to the LEFT the peak(s).')
            [X1,Y1] = ginput(1);
            figure(1);subplot(2,1,1);xlabel('Now click on the baseline to the RIGHT the peak(s).')
            [X2,Y2] = ginput(1);
            n=length(xx);
            %  Create "start" for this number of peaks
            yy=yy-((Y2-Y1)/(X2-X1)*(xx-X1)+Y1);
            [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
        case 8
            % When 'Backspace' key is pressed, user clicks the graph
            % along the presumed background points, then the program
            % subtracts the interploated background between the points.
            SavedSignal=Y;
            BaselinePoints=input('Number of baseline points to click): ');
            if isempty(BaselinePoints),BaselinePoints=8;end
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
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 92
            % When '\' key is pressed, restoreed original signal
            Y=SavedSignal;
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case {49,50,51,52,53,54,55,56,57}
            % When a number key is pressed, sets the number of peaks to that number.
            n=key-48;          
            ShapeString='';
            if round(n)~=NumPeaks,
                NumPeaks=round(n);
                FitResults=zeros(NumPeaks,5);
                ShapeString=SelectShapeString(Shape);
                subplot(2,1,2)
                if BIPOLAR,
                     xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = + - ' ShapeString  ] )
                else
                     xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = + ' ShapeString  ] )
                end               
            end % if
            n=max(xx)-min(xx);
            % Create a start value for this number of peaks
            start=[];
            startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker);
                start=[start markx n/(5*NumPeaks)];
            end
            PEAKHEIGHTS=zeros(NumPeaks,1);
            AA=zeros(NumPeaks,200);
            RedrawSignal(X,Y,xo,dx,NumPeaks);
            if BIPOLAR,
                xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = + - ' ShapeString  ] )
             else
                xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = + ' ShapeString  ] )
            end
        case 48  % Zero key is pressed to enter number of peaks
            disp(['Current number of peaks: ' num2str(NumPeaks)])
            NumPeaksInput=input('Type number or press Enter to keep unchanged: ');
            if isempty(NumPeaksInput),
            else
                NumPeaks=NumPeaksInput;
            end
            ShapeString=SelectShapeString(Shape);
            % Create a start value for this number of peaks
            n=max(xx)-min(xx);
            start=[];
            startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker);
                start=[start markx n/(5*NumPeaks)];
            end
            PEAKHEIGHTS=zeros(NumPeaks,1);
            AA=zeros(NumPeaks,200);
            RedrawSignal(X,Y,xo,dx,NumPeaks);
            if BIPOLAR,
                xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = + - ' ShapeString  ] )
            else
                xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = + ' ShapeString  ] )
            end
        case {103,108,111,112,101,104,59,106,117,115,71,76,96,72,66,80,123,69,85,86,68,84}
            % Selects peak shape when the following keys are pressed.
            switch key
                case 103 % When 'g' key is pressed, peak shape is set to Gaussian.
                    n=1;
                case 108 % When 'l' key is pressed, peak shape is set to Lorentzian.
                    n=2;
                case 111 % When 'o' key is pressed, peak shape is set to Lognormal distribution.
                    n=3;
                case 112 % When 'p' key is pressed, peak shape is set to Pearson.
                    n=4;
                case 101 % When 'e' key is pressed, peak shape is set to exponentally-broadened gaussian.
                    n=5;
                case 104 % When 'h' key is pressed, peak shape is set to equal width gaussians.
                    n=6;
                case 59 % When ';' key is pressed, peak shape is set to equal width gaussians.
                    n=7;
                case 106 % When 'J' key is pressed, peak shape is set to exponentally-broadened equal width gaussians.
                    n=8;
                case 117 % When 'U' key is pressed, peak shape is set to exponential pulse.
                    n=9;
                case 115 % When 'S' key is pressed, peak shape is set to Up Sigmoid (logistic function).
                    n=10;
                case 68  % When 'Shift-D' key is pressed, peak shape is set to Down Sigmoid (logistic function).
                    n=23;
                case 71 % When 'Shift-G' key is pressed, peak shape is set to FWGaussian.
                    n=11;
                    inputwidth=input('Peak width: ');
                    if isempty(inputwidth),
                    else
                        FIXEDPARAMETERS=inputwidth;
                    end
                case 76 % When ';' key is pressed, peak shape is set to FWLorentzian.
                    n=12;
                    inputwidth=input('Peak width: ');
                    if isempty(inputwidth),
                    else
                        FIXEDPARAMETERS=inputwidth;
                    end
                case 96 % When '`' key is pressed, peak shape is set to Gauss/Lorentz blend
                    n=13;
                case 72 % When 'Shift-H' key is pressed, peak shape is set to bifurcated Gaussian
                    n=14;
                case 66 % When Shift-B key is pressed, peak shape is set to Breit-Wigner-Fano
                    n=15;
                case 80 % When 'Shift-P' key is pressed, peak shape is set to Fixed-position Gaussians
                    n=16;
                    inputpositions=input('Peak positions as a vector, e.g. [200 400 600]: ');
                    if isempty(inputpositions),
                    else
                        FIXEDPARAMETERS=inputpositions;
                        FIXEDPOSITIONS=inputpositions;
                    end
                case 123 % When 'Shift-[' key is pressed, peak shape is set to Fixed-position Lorentzians
                    n=17;
                    inputpositions=input('Peak positions as a vector, e.g. [200 400 600]: ');
                    if isempty(inputpositions),
                    else
                        FIXEDPARAMETERS=inputpositions;
                    end
                case 69 % When 'Shift-E' key is pressed, peak shape is set to ExpLorentzian
                    n=18;
                case 85 % When 'Shift-U' key is pressed, peak shape is set to alpha function
                    n=19;
                case 86 % When 'Shift-V' key is pressed, peak shape is set to Voigt profile
                    n=20;   
                case 84  % When 'Shift-T' key is pressed, peak shape is set to triangular
                    n=21;
                otherwise
            end
            % switch
            if round(n)~=Shape,
                Shape=round(n);
                shapesvector=Shape;
                ShapeString=SelectShapeString(Shape);
                figure(1);subplot(2,1,2)
                if BIPOLAR,
                    xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = + - ' ShapeString  ] )
                else
                    xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = + ' ShapeString  ] )
                end
            end % if
        case 45 % if '-' key pressed
            disp(' ')
            disp('Select the peak shape of the model (type 1-20 and press Enter key):')
            disp('For a multi-shape model, type a vector of shape numbers,')
            disp('with brackets, e.g. [2 2 1 3]):')
            disp('Gaussians: y=exp(-((x-pos)./(0.6005615.*width)) .^2)')
            disp('  Gaussians with independent positions and widths : 1 (default)')
            disp('  Exponentional-broadened Gaussian (equal time constants): 5 ')
            disp('  Exponentional-broadened equal-width Gaussian : 8')
            disp('  Fixed-width exponentionally-broadened Gaussian = 36')
            disp('  Exponentional-broadened Gaussian (independent time constants): 31 ')
            disp('  Gaussians with the same widths : 6')
            disp('  Gaussians with preset fixed widths : 11')
            disp('  Fixed-position Gaussians : 16 ')
            disp('  Flattened Gaussian : 42 ')
            disp('  Asymmetrical Gaussians with unequal half-widths on both sides : 14')      
            disp('Lorentzians: y=ones(size(x))./(1+((x-pos)./(0.5.*width)).^2)')
            disp('  Lorentzians with independent positions and widths : 2')
            disp('  Exponentional-broadened Lorentzian : 18 ')            
            disp('  Equal-width Lorentzians : 7')
            disp('  Fixed-width Lorentzian : 12')
            disp('  Fixed-position Lorentzian : 17')
            disp('  Asymmetrical Lorentzians with unequal half-widths on both sides : 15')
            disp('Gaussian/Lorentzian blend (equal blends): 13')
            disp('  Fixed-width Gaussian/Lorentzian blend = 35')
            disp('  Gaussian/Lorentzian blend (independent blends): 33')
            disp('Voigt profile (fixed alphas): 20')
            disp('  Fixed-width Voigt profile (fixed alphas) : 34')
            disp('  Voigt profile (independent alphas): 30')  
            disp('Logistic: n=exp(-((x-pos)/(.477.*wid)).^2); y=(2.*n)./(1+n) : 3  ')
            disp('Pearson: y=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m : 4')
            disp('  Fixed-width Pearson = 37')
            disp('  Pearson with independent shape factors, m : 32')            
            disp('Exponential pUlse: y=exp(-tau1.*x).*(1-exp(-tau2.*x)): 9')
            disp('  equal-width exponential pulses: 48')
            disp('  Alpha function: y=(x-spoint)./pos.*exp(1-(x-spoint)./pos); : 19')           
            disp('Sigmoidal shapes')
            disp('  Up Sigmoid (logistic function): y=.5+.5*erf((x-tau1)/sqrt(2*tau2)) : 10')
            disp('  Down Sigmoid y=.5-.5*erf((x-tau1)/sqrt(2*tau2) ): 23')
            disp('  3-parameter Gompertz: 43')
            disp('  4-parameter logistic: 45')
            disp('Triangular: 21')
            disp('Sine wave: 40')
            disp('Rectangular pulse: 41')
            disp('Exponential decay -exp(-k*t): 44')
            disp('Blackbody radiation: 47')
            disp('Baseline shapes:')
            disp('  Linear slope: 26')
            disp('  Quadratic Baseline: 46'); 
            disp('For multi-shape, type numbers in brackets: e.g. [1 1 2 26] ')
            disp(' ')
            Shapeinput=input('Peak shape number or [vector]: ');
            if isempty(Shapeinput),
            else
                Shape=Shapeinput;
            end
            if isscalar(Shape),
            else
                % disp('peakshape is vector');
                shapesvector=Shape;
                NumPeaks=length(shapesvector);
                Shape=22;
                n=max(xx)-min(xx);
                extra=n/100.*zeros(size(shapesvector));
                start=[];
                startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
                for marker=1:NumPeaks,
                    markx=startpos(marker);
                    start=[start markx n/(5*NumPeaks)];
                end
                % start=start
            end
            if Shape>37, Shape=37;end
            if Shape<1, Shape=1;end
            if isempty(Shape),Shape=1;end
            switch Shape
                case 1
                    ShapeString='Gaussian';
                case 2
                    ShapeString='Lorentzian';
                case 3
                    ShapeString='Lognormal distribution';
                case 4
                    ShapeString='Pearson';
                case 5
                    ShapeString='ExpGaussian';
                case 6
                    ShapeString='Equal-width Gaussian';
                case 7
                    ShapeString='Equal-width Lorentzian';
                case 8
                    ShapeString='Equal-width ExpGauss.';
                case 9
                    ShapeString='Exponental pulse';
                case 10
                    ShapeString='Up Sigmoid (logistic function)';
                case 11
                    ShapeString='Fixed-width Gaussian';
                    inputwidth=input('Peak width vector [width1 width2...]: ');
                    if isempty(inputwidth),
                    else
                        FIXEDPARAMETERS=inputwidth;
                    end
                case 12
                    ShapeString='Fixed-width Lorentzian';
                    inputwidth=input('Peak width vector [width1 width2...]: ');
                    if isempty(inputwidth),
                    else
                        FIXEDPARAMETERS=inputwidth;
                    end
                case 13
                    ShapeString='Gauss/Lorentz blend';
                case 14
                    ShapeString='bifurcated Gaussian';
                case 15
                    ShapeString='Breit-Wigner-Fano';
                case 16
                    ShapeString='Fixed-position Gaussians';
                    inputpositions=input('Peak positions as a vector, e.g. [200 400 600]: ');
                    if isempty(inputpositions),
                    else
                        FIXEDPARAMETERS=inputpositions;
                    end
                case 17
                    ShapeString='Fixed-position Lorentzians';
                    inputpositions=input('Peak positions as a vector, e.g. [200 400 600]: ');
                    if isempty(inputpositions),
                    else
                        FIXEDPARAMETERS=inputpositions;
                    end
                case 18
                    ShapeString='ExpLorentzian';
                case 19
                    ShapeString='Alpha function';
                case 20
                    ShapeString='Voigt profile';
                case 21
                    ShapeString='triangular';
                case 22
                    ShapeString=num2str(shapesvector);
                case 23
                    ShapeString='Down Sigmoid (logistic function)';
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
                case 34
                    ShapeString='Fixed-width Voigt';
                    inputwidth=input('Peak width vector [width1 width2...]: ');
                    if isempty(inputwidth),
                    else
                        FIXEDPARAMETERS=inputwidth;
                    end
                case 35
                    ShapeString='Fixed-width G/L blend';
                    inputwidth=input('Peak width vector [width1 width2...]: ');
                    if isempty(inputwidth),
                    else
                        FIXEDPARAMETERS=inputwidth;
                    end
                case 36;
                    ShapeString='Fixed-width ExpGaussian';
                    inputwidth=input('Peak width vector [width1 width2...]: ');
                    if isempty(inputwidth),
                    else
                        FIXEDPARAMETERS=inputwidth;
                    end
                case 37
                    ShapeString='Fixed-width Pearson';
                    inputwidth=input('Peak width vector [width1 width2...]: ');
                    if isempty(inputwidth),
                    else
                        FIXEDPARAMETERS=inputwidth;
                    end
                case 38
                    ShapeString='ExpLorentzian (var. time constant)'; 
                case 40
                    ShapeString='Sine wave';
                case 41
                    ShapeString='Rectangular pulse';
                case 42
                    ShapeString='Flattened Gaussian';
                case 43
                    ShapeString='3-parameter Gompertz.';
                case 44
                    ShapeString='1-exp(-k*t)';
                case 45
                    ShapeString='4-parameter logistic';
                case 46
                    ShapeString='Quadratic Baseline';
                case 47
                    ShapeString='Blackbody';
                case 48
                    ShapeString='Equal-width exp. pulse';
                otherwise
            end % switch
            disp([ ShapeString  ' shape selected'])
            Shape=Shape
            figure(1);subplot(2,1,2)
            if BIPOLAR,
                xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = + - ' ShapeString  ] )
            else
                xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = + ' ShapeString  ] )
            end
        case 88 % Shift-X
            disp(['Current value of ''extra'' parameter: ' num2str(extra) ] )
            ExtraInput=input('Type value, [vector] or press Enter to keep unchanged: ');
            if isempty(ExtraInput),
            else
                extra=ExtraInput;
            end
            if Shape==4||5||13||14||15||18, % Only do a re-fit if shape uses Extra parameter.
                [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
            else
                [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
            end
        case 97
            % When 'a' key is pressed, increases "extra" by 5%
            extra=extra+.05*extra;
            if extra==0, extra=.01;end
            if Shape==4||5||13||14||15||18, % Only do a re-fit if shape uses Extra parameter.
                [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
            else
                [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
            end
        case 122
            % When 'z' key is pressed, decreases "extra" by 5%
            extra=extra-.05*extra;
            if extra==0, extra=.01;end
            if Shape==4||5||13||14||15||18, % Only do a re-fit if shape uses Extra parameter.
                [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
            else
                [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
            end
        case 65
            % When 'A' (Shift-A) key is pressed, increases "extra" by 0.5%
            extra=extra+.005*extra;
            if extra==0, extra=.001;end
            if Shape==4||5||13||14||15||18, % Only do a re-fit if shape uses Extra parameter.
                [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
            else
                [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
            end
        case 90
            % When 'Z' (Shift-Z) key is pressed, decreases "extra" by 0.5%
            extra=extra-.005*extra;
            if extra==0, extra=.001;end
            if Shape==4||5||13||14||15||18, % Only do a re-fit if shape uses Extra parameter.
                [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
            else
                [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
            end  
        case 114
            % When 'r' key is pressed, prints out table of fit results
            disp(['Fitting Error = ' num2str(MeanFitError(1)) '%' ])
            if Shape==9||Shape==10||Shape==19, % Pulse, alpha, and Sigmoid only
                disp('         Peak#     Tau1         Height       Tau2          Area');
            else
                if Shape==0, % For future use
                    disp('         Peak#      Position       Height      Area');
                else
                    disp('         Peak#   Position       Height       Width         Area');
                end
            end
            disp(FitResults)
            if AUTOZERO==3,
                disp([ 'Baseline= ' num2str(PEAKHEIGHTS(1)) ]);
            end
        case 100
            % When 'd' key is pressed, prints out table of model peak data
            disp('     x          y')
            disp([xx' yy'])
            if AUTOZERO==3,
                bl=1;
            else
                bl=0;
            end
            switch NumPeaks,
                case 1
                    disp('     x           peak 1')
                    disp([xxx' PEAKHEIGHTS(1+bl)*AA(1,:)' ])
                case 2
                    disp('     x           peak 1        peak 2')
                    disp([xxx' PEAKHEIGHTS(1+bl)*AA(1,:)' PEAKHEIGHTS(2+bl)*AA(2,:)'])
                case 3
                    disp('     x           peak 1        peak 2       peak 3')
                    disp([xxx' PEAKHEIGHTS(1+bl)*AA(1,:)' PEAKHEIGHTS(2+bl)*AA(2,:)' PEAKHEIGHTS(3+bl)*AA(3,:)'])
                case 4
                    disp('     x           peak 1        peak 2       peak 3       peak 4')
                    disp([xxx' PEAKHEIGHTS(1+bl)*AA(1,:)' PEAKHEIGHTS(2+bl)*AA(2,:)' PEAKHEIGHTS(3+bl)*AA(3,:)'  PEAKHEIGHTS(4+bl)*AA(4,:)'])
                case 5
                    disp('     x           peak 1        peak 2       peak 3       peak 4       peak 5')
                    disp([xxx' PEAKHEIGHTS(1+bl)*AA(1,:)' PEAKHEIGHTS(2+bl)*AA(2,:)' PEAKHEIGHTS(3+bl)*AA(3,:)'  PEAKHEIGHTS(4+bl)*AA(4,:)'  PEAKHEIGHTS(5+bl)*AA(5,:)'])
                case 6
                    disp('     x           peak 1        peak 2       peak 3       peak 4       peak 5       peak 6')
                    disp([xxx' PEAKHEIGHTS(1+bl)*AA(1,:)' PEAKHEIGHTS(2+bl)*AA(2,:)' PEAKHEIGHTS(3+bl)*AA(3,:)'  PEAKHEIGHTS(4+bl)*AA(4,:)'  PEAKHEIGHTS(5+bl)*AA(5,:)'  PEAKHEIGHTS(6+bl)*AA(6,:)'])
            end
        case 61 % When '=' key is pressed
            if BIPOLAR==0;BIPOLAR=1;else BIPOLAR=0;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);   
            subplot(2,1,1); if BIPOLAR==0;ylabel('+ mode');else ylabel('+ - mode');end
        case 78 % When Shift-N key pressed, negates signal
            Y=-Y;
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);          
        case 116
            % When 't' key is pressed, steps through AUTOZERO modes
            AUTOZERO=AUTOZERO+1;
            if AUTOZERO==4,AUTOZERO=0;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 89
            
        case 107
            % When 'k' key is pressed, prints out table of keyboard commands
            disp('KEYBOARD CONTROLS (Version 12):')
            disp(' Pan signal left and right...Coarse: < and >')
            disp('                             Fine: left and right cursor arrow keys')
            disp('                             Nudge: [ ] ')
            disp(' Zoom in and out.............Coarse zoom: / and "  ')
            disp('                             Fine zoom: up and down cursor arrow keys')
            disp(' Resets pan and zoom.........ESC')
            disp(' Select entire signal........Crtl-A')
            disp(' Select # of peaks...........Number keys 1-9, or press 0 key to enter manually')
            disp(' Select peak shape from menu - (minus or hyphen), type number or [vector] and Enter')
            disp(' Select peak shape by key....g Gaussian')
            disp('                             h equal-width Gaussians')
            disp('                             Shift-G  fixed-width Gaussians')
            disp('                             Shift P fixed-position Gaussians')
            disp('                             Shift-H  bifurcated Gaussians (a,z keys adjust shape)')
            disp('                             e Exponential-broadened Gaussian')
            disp('                             j exponential-broadened equal-width Gaussians')
            disp('                                 (a,z keys adjust broadening)')
            disp('                             l Lorentzian')
            disp('                             ; equal-width Lorentzians')
            disp('                             Shift [ fixed-position Lorentzians')
            disp('                             Shift-B Breit-Wigner-Fano (a,z keys adjust Fano factor)')
            disp('                             Shift-E Exponential-broadened Lorentzians')
            disp('                             Shift-L Fixed-width Lorentzians')
            disp('                             o Lognormal distribution')
            disp('                             p Pearson (a,z keys adjust shape)')
            disp('                             u exponential pUlse  y=exp(-tau1.*x).*(1-exp(-tau2.*x))')
            disp('                             Shift-U Alpha function: y=(x-tau2)./tau1.*exp(1-(x-tau2)./tau1)')
            disp('                             s Up Sigmoid (logistic function): y=.5+.5*erf((x-tau1)/sqrt(2*tau2))')
            disp('                             Shift-D Down Sigmoid y=.5-.5*erf((x-tau1)/sqrt(2*tau2))')
            disp('                             ~ Gauss/Lorentz blend (a/z keys adjust shape)')
            disp('                             V (Shift-V) Voigt profile (a/z adjusts shape)')
            disp('                             Shift-T  Triangular')
            disp(' Fit.........................f Performs one fit of the specified model')
            disp(' Select autozero mode........t  selects none, linear, quadratic, or flat baseline mode')
            disp(' Monopolar/bipolar mode......= +  Flips between + peaks only and +/- peak mode')
            disp(' Toggle log y mode OFF/ON....m  Plot linear or log Y axis on lower graph')
            disp(' 2-point Baseline............b, then click left and right baseline')
            disp(' Set manual baseline.........Backspace, then click baseline at multiple points')
            disp(' Restore original baseline...\  to cancel previous background subtraction')
            disp(' Click start positions.......c, click on peak position for each peak')
            disp(' Type in start vector........C (Shift-C) Type or Paste start vector [p1 w1 p2 w2 ...]')
            disp(' Enter value of ''extra''....Shift-x, type value or [vector], press Enter.')
            disp(' Adjust ''extra'' up/down....a,z: 5% change; upper case A,Z: 0.5% change.')
            disp(' Print parameters & results..q')
            disp(' Print fit results only......r')
            disp(' Compute bootstrap stats.....v  Estimates standard deViations of parameters.')
            disp(' Test effect of Noise........n  Test effect of Noise by fitting subset of points')
            disp(' Plot signal in figure 2.....y')
            disp(' Print model data table......d')
            disp(' Refine fit..................x  Takes best of 10 trial fits (change in line 177)')
            disp(' Print peakfit function......w  Print equivalent peakfit function with all input arguments')
            disp(' Save Figure as png file.....Shift-S  Saves as Figure1.png, Figure2.png, etc.')
            disp(' Display current settings....Shift-? displays table of current values')
            disp(' Fit polynomial to segment...Shift-o  Asks for polynomial order')
            disp(' Enter minimun width.........Shift-W')
            disp(' Enter saturation maximum....Shift-M')
            disp(' Switch to iPeak.............Shift-Ctrl-P  transfer current signal to iPeak.m')
            disp(' Switch to iSignal...........Shift-Ctrl-S  transfer current signal to iSignal.m')
        case 120 % When 'x' key is pressed, calls peakfit to take best of 'NumTrials' trial fits
            NumTrials=10;
            center=(max(xx)+min(xx))/2;
            window=max(xx)-min(xx);
%            disp(['Calculating best of ' num2str(NumTrials) ' trial fits....' ])
             if Shape==16||Shape==17, % Fixed-position shapes
                    fixedparameters=FIXEDPARAMETERS;
                else
                     fixedparameters=FIXEDPARAMETERS;
             end 
            FirstGuess=[];
            for peaknumber=1:NumPeaks,
                FirstGuess=[FirstGuess FitResults(peaknumber,2) FitResults(peaknumber,4)];
            end
            tic
            % FirstGuess=FirstGuess          
            [FitResults,MeanFitError]=peakfit([X',Y'],center,window,NumPeaks,shapesvector,extra,NumTrials,FirstGuess,AUTOZERO,fixedparameters,1,BIPOLAR,MINWIDTH,delta);
            et=toc;
%             disp(['     Elapsed time = '  num2str(et) ' sec.   Percent Fitting Error ' num2str(MeanFitError)])
%             if Shape==9||Shape==10||Shape==19, % Pulse, alpha and Sigmoid only
%                 disp('         Peak#     Tau1         Height       Tau2          Area');
%             else
%                 if Shape==0, % For future use
%                     disp('         Peak#      Position       Height      Area');
%                 else
%                     disp('         Peak#   Position       Height       Width         Area');
%                 end
%             end
%             disp(FitResults)
%             figure(1)
            start=FirstGuess;
        case 119
            % When 'W' is pressed, prints out peakfit functions with arguments in command window
            center=(max(xx)+min(xx))/2;
            window=max(xx)-min(xx);
            FirstGuess=[];
            if size(FitResults)==[0,0],
                disp('Perform at least one fit before using this command')
            else
                for peaknumber=1:NumPeaks,
                    FirstGuess=[FirstGuess FitResults(peaknumber,2) FitResults(peaknumber,4)];
                end
                disp(' ')
                disp('Copy and Paste these functions, replacing "datamatrix" with your data variable:')
                disp(['ipf(datamatrix,' num2str(center) ',' num2str(window) ')']);
                disp(['[FitResults,FitError]=peakfit(datamatrix,' num2str(center) ',' num2str(window) ',' num2str(NumPeaks) ',' num2str(Shape) ',' num2str(extra) ',' num2str(NumTrials) ', [' num2str(FirstGuess) '], ' num2str(AUTOZERO) ', [' num2str(FIXEDPARAMETERS) '], 1 ,' num2str(BIPOLAR) ',' num2str(MINWIDTH) ',' num2str(delta) ',' num2str(CLIPHEIGHT) ')' ] );
            end
         case 113
            % When 'q' key is pressed, prints out fitting parameters
            ShapeString=SelectShapeString(Shape);
            disp('------------------------------------------------------------------')
            if Shape==22,
                AllShapes=[];
                for NumShape=1:length(shapesvector),
                    AllShapes=[AllShapes  SelectShapeString(shapesvector(NumShape)) ', ' ];
                end
                disp( ['Peak shapes = ' AllShapes] )
                disp(['Extra vector = ' num2str(extra)])
            else
                disp(['Peak Shape = ' ShapeString])
            end
            switch Shape
                case 4
                    disp(['Shape Constant = ' num2str(extra)])
                case {5,8,18}
                    disp(['Time Constant = ' num2str(extra)])
                case 13
                    disp(['Percent Gaussian = ' num2str(extra)])
                case 14
                    disp(['Asymmetry = ' num2str(extra)])
                case 15
                    disp(['Fano factor = ' num2str(extra)])
                case 20
                    disp(['Alpha = ' num2str(extra)])
            end
            if BIPOLAR,
                disp('Bipolar mode')
            else
                disp('Positive peaks only')
            end
            if logplot,
                 disp('Log mode')
            else
                 disp('Linear mode')
            end
            switch AUTOZERO,
                case 0
                    disp('No baseline correction')
                case 1
                    disp('Linear baseline correction')
                case 2
                    disp('Quadratic baseline correction')
                case 3
                    disp('Flat baseline mode')
            end
            disp(['Number of peaks = ' num2str(NumPeaks)])
            
            if Shape==11||Shape==12, disp(['Fixed Peak Width = ' num2str(FIXEDPARAMETERS)]), end
            disp(['Fitted x range = ' num2str(min(xx)) ' - ' num2str(max(xx)) ' (dx=' num2str(max(xx)-min(xx)) ')  (Center=' num2str((max(xx)+min(xx))/2) ')  ' ])
            apnt=1;
            for pnt=1:length(xx),
                if yy(pnt)<CLIPHEIGHT,
                    axx(apnt)=xx(pnt);
                    ayy(apnt)=yy(pnt);
                    apnt=apnt+1;
                end
            end
            xx=axx;yy=ayy;
            disp([num2str(length(xx)) ' data points fit' ])
            disp(['Saturation height = ' num2str(CLIPHEIGHT) ])
            disp(['Minimum width = ' num2str(MINWIDTH) ])
            disp(['Percent Fitting Error = ' num2str(MeanFitError) '% ' ])
            if Shape==9||Shape==10||Shape==19,
                disp('         Peak#     Tau1         Height       Tau2          Area');
            else
                disp('         Peak#   Position       Height       Width         Area');
            end
            disp(FitResults)
            if AUTOZERO==3,
                disp([ 'Baseline= ' num2str(PEAKHEIGHTS(1)) ]);
            end
        case 121    % When 'Y' kkey is pressed (Added on version 5)
            figure(2) % Plot the entire signal cleanly in Figure window 2
            plot(X,Y)
            axis([X(1) X(length(X)) min(residual) max(Y)]);
            hold on
            for m=1:NumPeaks,
                % Add the individual component peaks in green lines
                if AUTOZERO==3,
                    plot(xxx,PEAKHEIGHTS(m+1)*AA(m,:),'g')
                else
                    plot(xxx,PEAKHEIGHTS(m)*AA(m,:),'g')
                end
            end
            % Show residual in red
            plot(xx,residual,'r')
            hold off
            title('Blue=Original Data; Green=Computed Fit; Red=Residual (Fit-Data)')
            % figure(1)
        case 118 % 'v' key computes bootstrap statistics (Added in version 8)
            NumTrialsBoot=input('Number of fit trials per bootstrap sample (0 to cancel): ');
            if isempty(NumTrialsBoot),NumTrialsBoot=1;end
            % EstTime=2.4+NumTrialsBoot.*(length(xx)./1436).*NumPeaks;
            EstTime=round((0.71581.*NumTrialsBoot).*(0.00069642.*length(xx)).*(5.4659.*NumPeaks)+2.5);
            % lengthxx=length(xx)
            % disp(['Estimated time: ',num2str(EstTime) ' seconds']);
            if NumTrialsBoot,
                disp('Computing bootstrap sampling statistics....May take several minutes.')
                tic;
                BootstrapResultsMatrix=zeros(5,100,NumPeaks);
                BootstrapErrorMatrix=zeros(2,100,NumPeaks);
                center=(max(xx)+min(xx))/2;
                window=max(xx)-min(xx);
                clear bx by
                cutoff=0.5;
                for trial=1:100,
                    n=1;
                    bx=X';
                    by=Y';
                    while n<length(X)-1,
                        if rand>cutoff,
                            bx(n)=X(n+1);
                            by(n)=Y(n+1);
                        end
                        n=n+1;
                    end
                    [FitResults,BootFitError]=peakfit([bx,by],center,window,NumPeaks,Shape,extra,NumTrialsBoot,start,AUTOZERO,FIXEDPARAMETERS,1,BIPOLAR,MINWIDTH,delta);
                    for peak=1:NumPeaks,
                        BootstrapResultsMatrix(:,trial,peak)=FitResults(peak,:);
                        BootstrapErrorMatrix(:,trial,peak)=BootFitError;
                    end
                end
                for peak=1:NumPeaks,
                    disp(' ')
                    if Shape==9||Shape==10||Shape==19, % Pulse, alpha and Sigmoid only
                        disp(['     Peak #' num2str(peak) '    Tau1         Height       Tau2          Area']);
                    else
                        if Shape==0, % For future use
                            disp(['     Peak #' num2str(peak) '    Position       Height      Area']);
                        else
                            disp(['     Peak #' num2str(peak) '    Position       Height       Width         Area']);
                        end
                    end
                    BootstrapMean=mean(real(BootstrapResultsMatrix(:,:,peak)'));
                    BootstrapSTD=std(BootstrapResultsMatrix(:,:,peak)');
                    BootstrapIQR=iqr(BootstrapResultsMatrix(:,:,peak)');
                    PercentRSD=100.*BootstrapSTD./BootstrapMean;
                    PercentIQR=100.*BootstrapIQR./BootstrapMean;
                    MaxError=max(real(BootstrapErrorMatrix(:,:,peak)'));
                    MinError=min(real(BootstrapErrorMatrix(:,:,peak)'));
                    disp(['Mean:        ', num2str(BootstrapMean(2:5))])
                    disp(['STD:         ', num2str(BootstrapSTD(2:5))])
                    disp(['RSD (IQR):   ', num2str(BootstrapIQR(2:5))])
                    disp(['% RSD:       ', num2str(PercentRSD(2:5))])
                    disp(['% RSD (IQR): ', num2str(PercentIQR(2:5)) '% ' ])
                    % mean(PercentIQR(2:5)./PercentRSD(2:5))
                end
                toc;
                disp(['Max/Min Fitting Error of the 100 bootstrap samples: ', num2str(MaxError) '/' num2str(MinError) ])
                disp('-------------------------------------------------------------------')
            end
            figure(1)
            title('ipf 12.2  Typical Bootstrap sample fit')
        case 110 % When 'n' key is pressed (Added on version 8)
%            disp('Fit to single bootstrap sample')
            tic;
            center=(max(xx)+min(xx))/2;
            window=max(xx)-min(xx);
            clear bx by
            cutoff=0.5;
            n=1;
            bx=X';
            by=Y';
            while n<length(X)-1,
                if rand>cutoff,
                    bx(n)=X(n+1);
                    by(n)=Y(n+1);
                end
                n=n+1;
            end
            StartVector=start;
            startnow=start;
            delta=(max(xx)-min(xx))/100;
            for k=1:2*NumPeaks,
                startnow(k)=start(k)+(rand-.5)*delta;
            end
            start=startnow; 
            [FitResults,MeanFitError]=peakfit([bx,by],center,window,NumPeaks,shapesvector,extra,NumTrialsBoot,start,AUTOZERO,FIXEDPARAMETERS,1,BIPOLAR,MINWIDTH,delta);
%             disp(['Fitting Error = ' num2str(MeanFitError) '%'])
%             if Shape==9||Shape==10||Shape==19,
%                 disp('         Peak#     Tau1         Height       Tau2          Area');
%             else
%                 disp('         Peak#   Position       Height       Width         Area');
%             end
%             disp(FitResults)
            figure(1)
            title('ipf 12.2  Single Bootstrap sample fit')
        case 83
            saveas(gcf,['Figure' num2str(FigNum) '.png'])
            FigNum=FigNum+1;
        case 63 % Shift-? Displays current settings
            disp(' ')
            disp( 'Current Settings:')
            disp([ 'Fitted x range = ' num2str(min(xx)) ' - ' num2str(max(xx)) ' (dx=' num2str(max(xx)-min(xx)) ')  (Center=' num2str((max(xx)+min(xx))/2) ')  ' ])
            
            disp([ 'NumPeaks = ' num2str(NumPeaks) ])
            disp([ 'Shape = ' num2str(Shape) ])
            disp([ 'shapesvector = ' num2str(shapesvector) ])
            disp([ 'extra = ' num2str(extra) ])
            disp([ 'NumTrials = ' num2str(NumTrials) ])           
            disp([ 'start vector = [' num2str(start) ']' ])
            disp([ 'Baseline mode = ' num2str(AUTOZERO) ])
            disp([ 'Fixed parameters = ' num2str(FIXEDPARAMETERS) ])  
            disp([ 'BIPOLAR = ' num2str(BIPOLAR) ])
            disp([ 'Minimum peak width = ' num2str(MINWIDTH) ])
            disp([ 'delta = ' num2str(delta) ])            
            disp([ 'Saturation height = ' num2str(CLIPHEIGHT) ])

        case 19 % Shift-Ctrl-F transfers current signal to iSignal
            isignal(X,Y);
        case 16 % Shift-Ctrl-P transfers current signal to Interactive Peak Detector
            ipeak(X,Y);
        case 79 % Shift-o
            polyorder=input('Polynomial order (1=linear, 2=quadratic, etc): ');
            if isempty(polyorder),polyorder=0;end
            if polyorder>0,
                figure(1);subplot(2,1,1); % Select upper window
                [coef, RSquared]=plotit(xx,yy,polyorder);
                residual=yy-polyval(coef,xx);
                subplot(2,1,2);plot(xx,residual,'r.')
                xlabel('Residuals')
                disp(' ')
                disp([ 'Selected x range: ' num2str(length(xx)) ' points from x = ' num2str(min(xx)) ' to ' num2str(max(xx)) ])
                disp([ 'Polynomial order: ' num2str(polyorder)])
                disp([ 'Polynomial coefficients (slope, intercept): ' num2str(coef)])
                disp([ 'Coefficient of determination (R2): ' num2str(RSquared)])
            end
        case 87 % Set minimum peak width
            disp(['Current value of minimun width: ' num2str(MINWIDTH) ] )
            MinWidthInput=input('Type value or press Enter to keep unchanged: ');
            if isempty(MinWidthInput),
            else
                MINWIDTH=MinWidthInput;
            end   
            figure(1)
        case 77 % Shift-M sets Saturation height, ignores y values above this value
            disp(['Current value Saturation height: ' num2str(CLIPHEIGHT) ] )
            ClipInput=input('Type value or press Enter to keep unchanged: ');
            if isempty(ClipInput),
            else
                CLIPHEIGHT=ClipInput;
            end 
            figure(1)
        otherwise
            UnassignedKey=double(key)
            disp('Press k to print out list of keyboard commands')
            disp('Make sure the CAPS LOCK mode is not engaged.')
    end % switch
end % if
% ----------------------------------------------------------------------
function [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks)
% Plots the entire signal (X,Y) in the lower half of the plot window and an
% isolated segment (xx,yy) in the upper half, controlled by Pan and Zoom
% keys. Top half of the figure shows original signal
global AUTOZERO AA xxx PEAKHEIGHTS logplot BIPOLAR MINWIDTH
minX=min(X);maxX=max(X);minY=min(Y);maxY=max(Y);
Startx=round(xo-(dx/2));
Endx=abs(round(xo+(dx/2)-1));
if Endx>length(Y),Endx=length(Y);end
if Startx<1,Startx=1;end
PlotRange=Startx:Endx;
if (Endx-Startx)<2, PlotRange=xo:xo+2;end
xx=X(PlotRange);
yy=Y(PlotRange);
X1=min(xx);
X2=max(xx);
Y1=min(Y);
Y2=max(Y);
% Remove baseline from data segment
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

hold off
clf
figure(1);subplot(2,1,1); % Select upper window
% semilogy(xx,yy,'b.');
if logplot,
    semilogy(xx,yy,'b.');
else
    plot(xx,yy,'b.'); % Plot the original signal in blue in upper window
end
xlabel('Line up the estimated peak positions roughly with the vertical lines')
lyy=min(yy);
% lyy=min(yy);
uyy=max(yy)+(max(yy)-min(yy))/10;
switch AUTOZERO,
    case 0
        title('ipf 12.2  No baseline correction. Pan and Zoom to isolate peaks to be fit in upper window.')
        if lyy<uyy;axis([X(Startx) X(Endx) lyy uyy ]);end
    case 1
        title('ipf 12.2  Linear baseline mode. Pan and Zoom to isolate peaks to be fit in upper window.')
        if lyy<uyy;axis([X(Startx) X(Endx) lyy uyy ]);end
    case 2
        title('ipf 12.2  Quadratic baseline mode. Pan and Zoom to isolate peaks to be fit in upper window.')
        if lyy<uyy;axis([X(Startx) X(Endx) lyy uyy ]);end
    case 3
        title('ipf 12.2  Flat baseline mode. Pan and Zoom to isolate peaks to be fit in upper window.')
        if lyy<uyy;axis([X(Startx) X(Endx) lyy uyy ]);end
end
hold on
% Mark starting peak positions with vertical dashed lines in upper window
% Determine locations for peak (vertical line) markers
n=X2-X1;
width=n/(5*NumPeaks);
start=[];
startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+X1;
for marker=1:NumPeaks,
    markx=startpos(marker);
    start=[start markx width];
    semilogy([markx markx],[lyy uyy],'m--')
end % for marker
hold off
%
% Bottom half of the figure shows full signal in either linear or log mode
% as set by the M key.
subplot(2,1,2);cla
if logplot,
    semilogy(X,abs(Y))  % Graph the signal with linear Y axis
    ylabel('Log y mode')
    axis([X(1) X(length(X)) min(abs(Y)) max(Y)]); % Update plot
else
    plot(X,Y)  % Graph the signal with linear Y axis
    ylabel('Linear y mode')
    axis([X(1) X(length(X)) min(Y) max(Y)]); % Update plot
end
title('# peaks 1-9,0   Shapes: g h G H e j l ; L : o p u s `   Fit=f   t=autozero   m=linear/log')
hold on
for marker=1:NumPeaks,
    markx=startpos(marker);
    plot([markx markx],[minY maxY],'m--')
end % for marker

% Mark the limits of the upper windows on the lower whole-signal plot
plot([X1 X1],[minY maxY],'g--')
plot([X2 X2],[minY maxY],'g--')
hold off
xlabel(['Press K to print out all keyboard commands'])
% ----------------------------------------------------------------------
function [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra)
% Given isolated segment of signal (xx,yy), plots it in upper half, computes fit with
% "NumPeaks" component peaks of shape "Shape", starting with start values
% "start", then plots residuals in lower half.
%  T. C. O'Haver (toh@umd.edu),  Version 2.2,  October, 2011.
global PEAKHEIGHTS AUTOZERO AA xxx residual FIXEDPARAMETERS BIPOLAR logplot
global shapesvector MINWIDTH CLIPHEIGHT
PEAKHEIGHTS=zeros(1,NumPeaks);
% Saturation Code, skips points greater than set maximum
apnt=1;
for pnt=1:length(xx),
     if yy(pnt)<CLIPHEIGHT,
         axx(apnt)=xx(pnt);
         ayy(apnt)=yy(pnt);
         apnt=apnt+1;
     end
end
xx=axx;yy=ayy;
n=length(xx);
x0=min(xx);
% StartVector=start
% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.000001,'TolFun',.00001,'Display','off','MaxFunEvals',1000);
tic
switch Shape
    case 1
        FitParameters=fminsearch(@fitgaussian,start,options,xx,yy);
        ShapeString='Gaussian';
        for Peak=1:NumPeaks;
            if FitParameters(2*Peak)<MINWIDTH,
                FitParameters(2*Peak)=MINWIDTH;
            end
        end
    case 2
        FitParameters=fminsearch(@fitlorentzian,start,options,xx,yy);
        ShapeString='Lorentzian';
        for Peak=1:NumPeaks;
            if FitParameters(2*Peak)<MINWIDTH,
                FitParameters(2*Peak)=MINWIDTH;
            end
        end
    case 3
        FitParameters=fminsearch(@fitLognormal,start,options,xx,yy);
        ShapeString='Lognormal';
        for Peak=1:NumPeaks;
            if FitParameters(2*Peak)<MINWIDTH,
                FitParameters(2*Peak)=MINWIDTH;
            end
        end
    case 4
        FitParameters=fminsearch(@fitpearson,start,options,xx,yy,extra(1));
        ShapeString='Pearson';
        for Peak=1:NumPeaks;
            if FitParameters(2*Peak)<MINWIDTH,
                FitParameters(2*Peak)=MINWIDTH;
            end
        end
    case 5
        FitParameters=fminsearch(@fitexpgaussian,start,options,xx,yy,-extra(1));
        ShapeString='ExpGaussian';
        for Peak=1:NumPeaks;
            if FitParameters(2*Peak)<MINWIDTH,
                FitParameters(2*Peak)=MINWIDTH;
            end
        end
    case 6
        cwnewstart(1)=start(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=start(2.*pc-1);
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        FitParameters=fminsearch(@fitewgaussian,cwnewstart,options,xx,yy);
        ShapeString='Equal width Gauss.';
        for Peak=1:NumPeaks;
            if FitParameters(NumPeaks+1)<MINWIDTH,
                FitParameters(NumPeaks+1)=MINWIDTH;
            end
        end
    case 7
        cwnewstart(1)=start(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=start(2.*pc-1);
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        FitParameters=fminsearch(@fitewlorentzian,cwnewstart,options,xx,yy);
        ShapeString='Equal width Lorentzian';
        for Peak=1:NumPeaks;
            if FitParameters(NumPeaks+1)<MINWIDTH,
                FitParameters(NumPeaks+1)=MINWIDTH;
            end
        end
    case 8
        cwnewstart(1)=start(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=start(2.*pc-1);
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        FitParameters=fminsearch(@fitexpewgaussian,cwnewstart,options,xx,yy,-extra(1));
        ShapeString='Exp. equal width Gaussians';
        for Peak=1:NumPeaks;
            if FitParameters(NumPeaks+1)<MINWIDTH,
                FitParameters(NumPeaks+1)=MINWIDTH;
            end
        end
    case 9
        FitParameters=fminsearch(@fitexppulse,start,options,xx,yy);
        ShapeString='Exponential Pulse';
        for Peak=1:NumPeaks;
            if FitParameters(2*Peak)<MINWIDTH,
                FitParameters(2*Peak)=MINWIDTH;
            end
        end
    case 10
        FitParameters=fminsearch(@fitupsigmoid,start,options,xx,yy);
        ShapeString='Up Sigmoid (logistic function)';
        for Peak=1:NumPeaks;
            if FitParameters(2*Peak)<MINWIDTH,
                FitParameters(2*Peak)=MINWIDTH;
            end
        end
    case 23
        FitParameters=fminsearch(@fitdownsigmoid,start,options,xx,yy);
        ShapeString='Down Sigmoid (logistic function)';
        for Peak=1:NumPeaks;
            if FitParameters(2*Peak)<MINWIDTH,
                FitParameters(2*Peak)=MINWIDTH;
            end
        end
    case 11
        fixedstart=[];
        for pk=1:NumPeaks,
            fixedstart(pk)=start(2*pk-1);
        end
        % fixedstart=fixedstart % Testing
        FitParameters=fminsearch(@FitFWGaussian,fixedstart,options,xx,yy);
        ShapeString='Fixed-width Gaussian';
    case 12
        fixedstart=[];
        for pk=1:NumPeaks,
            fixedstart(pk)=start(2*pk-1);
        end
        FitParameters=fminsearch(@FitFWLorentzian,fixedstart,options,xx,yy);
        ShapeString='Fixed-width Lorentzian';
    case 13
        FitParameters=fminsearch(@fitGL,start,options,xx,yy,extra(1));
        ShapeString='Gausss/Lorentz blend';
    case 14
        FitParameters=fminsearch(@fitBiGaussian,start,options,xx,yy,extra(1));
        ShapeString='bifurcated Gaussian';
        for Peak=1:NumPeaks;
            if FitParameters(2*Peak)<MINWIDTH,
                FitParameters(2*Peak)=MINWIDTH;
            end
        end
    case 15
        FitParameters=fminsearch(@fitBWF,start,options,xx,yy,extra(1));
        ShapeString='Breit-Wigner-Fano';
        for Peak=1:NumPeaks;
            if FitParameters(2*Peak)<MINWIDTH,
                FitParameters(2*Peak)=MINWIDTH;
            end
        end
    case 16
        fixedstart=[];
        for pc=1:NumPeaks,
            fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
            fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc);
        end
        FitParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart,options);
        ShapeString='Fixed-position Gaussians';
        for Peak=1:NumPeaks;
            if FitParameters(Peak)<MINWIDTH,
                FitParameters(Peak)=MINWIDTH;
            end
        end
    case 17
        fixedstart=[];
        for pc=1:NumPeaks,
            fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
            fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc);
        end
        FitParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart,options);
        ShapeString='Fixed-position Lorentzians';
        for Peak=1:NumPeaks;
            if FitParameters(Peak)<MINWIDTH,
                FitParameters(Peak)=MINWIDTH;
            end
        end
    case 18
        FitParameters=fminsearch(@fitexplorentzian,start,options,xx,yy,-extra(1));
        ShapeString='ExpLorentzian';
        for Peak=1:NumPeaks;
            if FitParameters(Peak)<MINWIDTH,
                FitParameters(Peak)=MINWIDTH;
            end
        end
    case 19
        FitParameters=fminsearch(@fitalphafunction,start,options,xx,yy);
        ShapeString='Alpha function';
        for Peak=1:NumPeaks;
            if FitParameters(Peak)<MINWIDTH,
                FitParameters(Peak)=MINWIDTH;
            end
        end
    case 20
        FitParameters=fminsearch(@fitvoigt,start,options,xx,yy,extra(1));
        ShapeString='Voigt profile';
        for Peak=1:NumPeaks;
            if FitParameters(Peak)<MINWIDTH,
                FitParameters(Peak)=MINWIDTH;
            end
        end
    case 21
        FitParameters=fminsearch(@(lambda)(fittriangular(lambda,xx,yy)),start,options);
        ShapeString='Triangular';
        for Peak=1:NumPeaks;
            if FitParameters(Peak)<MINWIDTH,
                FitParameters(Peak)=MINWIDTH;
            end
        end
    case 22
        FitParameters=fminsearch(@(lambda)(fitmultiple(lambda,xx,yy,shapesvector,extra)),start,options);
        ShapeString=num2str(shapesvector);
        for Peak=1:NumPeaks;
            if FitParameters(Peak)<MINWIDTH,
                FitParameters(Peak)=MINWIDTH;
            end
        end
    otherwise
end % switch
et=toc;
% Construct model from fitted parameters
A=zeros(NumPeaks,n);
AA=zeros(NumPeaks,200);
xxx=linspace(min(xx),max(xx),200);
for m=1:NumPeaks,
    switch Shape
        case 1
            A(m,:)=gaussian(xx,FitParameters(2*m-1),FitParameters(2*m));
            AA(m,:)=gaussian(xxx,FitParameters(2*m-1),FitParameters(2*m));
        case 2
            A(m,:)=lorentzian(xx,FitParameters(2*m-1),FitParameters(2*m));
            AA(m,:)=lorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m));
        case 3
            A(m,:)=Lognormal(xx,FitParameters(2*m-1),FitParameters(2*m));
            AA(m,:)=Lognormal(xxx,FitParameters(2*m-1),FitParameters(2*m));
        case 4
            A(m,:)=pearson(xx,FitParameters(2*m-1),FitParameters(2*m),extra(1));
            AA(m,:)=pearson(xxx,FitParameters(2*m-1),FitParameters(2*m),extra(1));
        case 5
            A(m,:)=expgaussian(xx,FitParameters(2*m-1),FitParameters(2*m),-extra(1))';
            AA(m,:)=expgaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra(1)*length(xxx)./length(xx))';
        case 6
            A(m,:)=gaussian(xx,FitParameters(m),FitParameters(NumPeaks+1));
            AA(m,:)=gaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
        case 7
            A(m,:)=lorentzian(xx,FitParameters(m),FitParameters(NumPeaks+1));
            AA(m,:)=lorentzian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
        case 8
            A(m,:)=expgaussian(xx,FitParameters(m),FitParameters(NumPeaks+1),-extra(1));
            AA(m,:)=expgaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1),-extra(1)*length(xxx)./length(xx)');
        case 9
            A(m,:)=exppulse(xx,FitParameters(2*m-1),FitParameters(2*m));
            AA(m,:)=exppulse(xxx,FitParameters(2*m-1),FitParameters(2*m));
        case 10
            A(m,:)=upsigmoid(xx,FitParameters(2*m-1),FitParameters(2*m));
            AA(m,:)=upsigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m));
        case 23
            A(m,:)=downsigmoid(xx,FitParameters(2*m-1),FitParameters(2*m));
            AA(m,:)=downsigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m));
        case 11
            A(m,:)=gaussian(xx,FitParameters(m),FIXEDPARAMETERS);
            AA(m,:)=gaussian(xxx,FitParameters(m),FIXEDPARAMETERS);
        case 12
            A(m,:)=lorentzian(xx,FitParameters(m),FIXEDPARAMETERS);
            AA(m,:)=lorentzian(xxx,FitParameters(m),FIXEDPARAMETERS);
        case 13
            A(m,:)=GL(xx,FitParameters(2*m-1),FitParameters(2*m),extra(1));
            AA(m,:)=GL(xxx,FitParameters(2*m-1),FitParameters(2*m),extra(1));
        case 14
            A(m,:)=BiGaussian(xx,FitParameters(2*m-1),FitParameters(2*m),extra(1));
            AA(m,:)=BiGaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra(1));
        case 15
            A(m,:)=BWF(xx,FitParameters(2*m-1),FitParameters(2*m),extra(1));
            AA(m,:)=BWF(xxx,FitParameters(2*m-1),FitParameters(2*m),extra(1));
        case 16
            A(m,:)=gaussian(xx,FIXEDPARAMETERS(m),FitParameters(m));
            AA(m,:)=gaussian(xxx,FIXEDPARAMETERS(m),FitParameters(m));
        case 17
            A(m,:)=lorentzian(xx,FIXEDPARAMETERS(m),FitParameters(m));
            AA(m,:)=lorentzian(xxx,FIXEDPARAMETERS(m),FitParameters(m));
        case 18
            A(m,:)=explorentzian(xx,FitParameters(2*m-1),FitParameters(2*m),-extra(1))';
            AA(m,:)=explorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra(1)*length(xxx)./length(xx))';
        case 19
            A(m,:)=alphafunction(xx,FitParameters(2*m-1),FitParameters(2*m));
            AA(m,:)=alphafunction(xxx,FitParameters(2*m-1),FitParameters(2*m));
        case 20
            A(m,:)=voigt(xx,FitParameters(2*m-1),FitParameters(2*m),extra(1));
            AA(m,:)=voigt(xxx,FitParameters(2*m-1),FitParameters(2*m),extra(1));
        case 21
            A(m,:)=triangular(xx,FitParameters(2*m-1),FitParameters(2*m));
            AA(m,:)=triangular(xxx,FitParameters(2*m-1),FitParameters(2*m));
        case 22
            A(m,:)=peakfunction(shapesvector(m),xx,FitParameters(2*m-1),FitParameters(2*m),extra(m));
            AA(m,:)=peakfunction(shapesvector(m),xxx,FitParameters(2*m-1),FitParameters(2*m),extra(m));
        otherwise
    end % switch
end % for
% PEAKHEIGHTS=PEAKHEIGHTS % Error testing
% SizeA=size(A) % Error testing
% SizeAA=size(AA) % Error testing
% model=PEAKHEIGHTS'*A;  % Multiplies each row by the corresponding amplitude and adds them up
% mmodel=PEAKHEIGHTS'*AA;
if AUTOZERO==3,
    baseline=PEAKHEIGHTS(1);
    Heights=PEAKHEIGHTS(2:1+NumPeaks);
    model=Heights'*A+baseline;
    mmodel=Heights'*AA+baseline;
else
    model=PEAKHEIGHTS'*A;
    mmodel=PEAKHEIGHTS'*AA;
    Heights=PEAKHEIGHTS;
    baseline=0;
end
% Top half of the figure shows original signal and the fitted model.
figure(1);subplot(2,1,1);
if logplot,
    semilogy(xx,yy,'b.'); % Plot the original signal in blue dots
else
    plot(xx,yy,'b.'); % Plot the original signal in blue dots
end
% lyy=0;
lyy=min(yy);
uyy=max(yy)+(max(yy)-min(yy))/10;
axis([min(xx) max(xx) lyy uyy]);
if BIPOLAR,
    ylabel('+ - mode')
else
    ylabel('+ mode')
end
hold on
for m=1:NumPeaks,
    if logplot,
        semilogy(xxx,Heights(m)*AA(m,:),'g')  % Plot the individual component peaks in green lines
    else
        plot(xxx,Heights(m)*AA(m,:),'g')  % Plot the individual component peaks in green lines
    end
    area(m)=trapz(xx,Heights(m)*A(m,:)); % Compute the area of each component peak using trapezoidal method
end
axis([min(xx) max(xx) lyy uyy]);
% Mark starting peak positions with vertical dashed lines
for marker=1:NumPeaks,
    markx=start((2*marker)-1);
    subplot(2,1,1);plot([markx markx],[0 max(yy)],'m--')
end % for
if logplot,
    semilogy(xxx,mmodel,'r');  % Plot the total model (sum of component peaks) in red lines
else
    plot(xxx,mmodel,'r');  % Plot the total model (sum of component peaks) in red lines
end
hold off;
% lyy=min(yy);
% uyy=max(yy)+(max(yy)-min(yy))/10;
% if BIPOLAR,
%     axis([min(xx) max(xx) lyy uyy]);
%     ylabel('+ - mode')
% else
%     axis([min(xx) max(xx) 0 uyy]);
%     ylabel('+ mode')
% end
switch AUTOZERO,
    case 0
        title('ipf 12.2 No baseline correction. Pan and Zoom to isolate peaks to be fit in upper window.')
    case 1
        title('ipf 12.2 Linear baseline mode. Pan and Zoom to isolate peaks to be fit in upper window.')
    case 2
        title('ipf 12.2 Quadratic baseline mode. Pan and Zoom to isolate peaks to be fit in upper window.')
    case 3
        title('ipf 12.2 Flat baseline mode. Pan and Zoom to isolate peaks to be fit in upper window.')
end
xlabel('Vertical dotted lines indicate first guess peak positions. C to customize.');
% Bottom half of the figure shows the residuals and displays RMS error
% between original signal and model
if logplot,
    residual=log10(yy)-log10(model);
else
    residual=yy-model;
end
figure(1);subplot(2,1,2);plot(xx,residual,'r.')
if logplot,
     MeanFitError=100*norm(residual)./(sqrt(n)*(log10(max(yy)-log10(min(yy)))));
else
     MeanFitError=100*norm(residual)./(sqrt(n)*max(yy));
end
title('F: fit, X: best of 10 trials, T: baseline, Shift-X: enter Extra, a/z adjust Extra, Q: report')
ylabel('Residuals')
switch Shape
    case 4
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) ' Shape Constant = ' num2str(extra) '   Error = ' num2str(round(10000*MeanFitError)/10000) '%' ] )
    case {5,8,18}
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '    Time Constant = ' num2str(extra(1)) '    Min. Width = ' num2str(MINWIDTH) '  Error = ' num2str(round(10000*MeanFitError)/10000) '%' ] )
    case 13
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '       % Gaussian = ' num2str(extra(1)) '    Min. Width = ' num2str(MINWIDTH) '  Error = ' num2str(round(10000*MeanFitError)/10000) '%' ] )
    case 14
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '        Asymmetry = ' num2str(extra(1)) '    Min. Width = ' num2str(MINWIDTH) '  Error = ' num2str(round(10000*MeanFitError)/10000) '%' ] )
    case 15
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '        Fano factor = ' num2str(extra(1)) '    Min. Width = ' num2str(MINWIDTH) '  Error = ' num2str(round(10000*MeanFitError)/10000) '%' ] )
    case 20
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '            Alpha = ' num2str(extra(1)) '     Min. Width = ' num2str(MINWIDTH) ' Error = ' num2str(round(10000*MeanFitError)/10000) '%' ] )
    otherwise
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      Error = ' num2str(round(10000*MeanFitError)/10000) '% ' ] )
end
minres=min(residual);
maxres=max(residual);
if minres<maxres,
    axis([min(xx) max(xx) minres maxres ]);
else
    yysize=size(yy);
    modelsize=size(model);
end

% Put results into a matrix, one row for each peak, showing peak index number,
% position, amplitude, and width.
FitResults=[];
for m=1:NumPeaks,
    if m==1,
        if Shape==6||Shape==7||Shape==8, % Equal-width shapes
            FitResults=[[round(m) FitParameters(m) Heights(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if Shape==11||Shape==12, % Fixed-width shapes
                FitResults=[[round(m) FitParameters(m) Heights(m) FIXEDPARAMETERS area(m)]];
            else
                if Shape==16||Shape==17, % Fixed-position shapes
                     FitResults=[[round(m) FIXEDPARAMETERS(m) Heights(m) FitParameters(m) area(m)]];
                else
                    FitResults=[[round(m) FitParameters(2*m-1) Heights(m) abs(FitParameters(2*m)) area(m)]];
                end
            end
        end % if shape
    else
        if Shape==6||Shape==7||Shape==8,  % Equal-width shapes
            FitResults=[FitResults ; [round(m) FitParameters(m) Heights(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if Shape==11||Shape==12, % Fixed-width shapes
                FitResults=[FitResults ; [round(m) FitParameters(m) Heights(m) FIXEDPARAMETERS area(m)]];
            else
                if Shape==16||Shape==17, % Fixed-position shapes
                     FitResults=[FitResults ; [round(m) FIXEDPARAMETERS(m) Heights(m) FitParameters(m) area(m)]]; 
                else
                    FitResults=[FitResults ; [round(m) FitParameters(2*m-1) Heights(m) abs(FitParameters(2*m)) area(m)]];
                end
            end
        end % if shape
    end % m==1
end % for m=1:NumPeaks

% Display Fit Results on lower  graph
subplot(2,1,2);
startx=min(xx)+(max(xx)-min(xx))./20;
dxx=(max(xx)-min(xx))./10;
dyy=((max(residual)-min(residual))./10);
starty=max(residual)-dyy;
FigureSize=get(gcf,'Position');
if Shape(1)==9||Shape(1)==10||Shape(1)==19,  % Pulse and sigmoid shapes only
    text(startx,starty+dyy/2,['Peak #          tau1           Height           tau2             Area'] );
else
    text(startx,starty+dyy/2,['Peak #          Position        Height         Width             Area'] );
end
% Display FitResults using sprintf
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
end
% ----------------------------------------------------------------------
function ShapeString=SelectShapeString(Shape)
global shapesvector
switch Shape
    case 1
        ShapeString='Gaussian';
    case 2
        ShapeString='Lorentzian';
    case 3
        ShapeString='Lognormal distribution';
    case 4
        ShapeString='Pearson7';
    case 5
        ShapeString='Exp. Gaussian';
    case 6
        ShapeString='Equal width Gaussians';
    case 7
        ShapeString='Equal width Lorentzians';
    case 8
        ShapeString='Equal-width ExpGauss.';
    case 9
        ShapeString='Exponental pulse';
    case 10
        ShapeString='Up Sigmoid (logistic function)';
    case 11
        ShapeString='Fixed-width Gaussian';
    case 12
        ShapeString='Fixed-width Lorentzian';
    case 13
        ShapeString='Gauss/Lorentz blend';
    case 14
        ShapeString='bifurcated Gaussian';
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
        ShapeString='Voigt profile';
    case 21
        ShapeString='triangular';
    case 22
        ShapeString=num2str(shapesvector);
    case 23
        ShapeString='Down Sigmoid (logistic function)';
    otherwise
        ShapeString='';
end % switch Shape
% ----------------------------------------------------------------------
function [FitResults,GOF,baseline,coeff,residual,xi,yi,BootResults]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots,bipolar,minwidth,DELTA,clipheight)
% A command-line peak fitting program for time-series signals, written as a
% self-contained Matlab function in a single m-file. Uses a non-linear
% optimization algorithm to decompose a complex, overlapping-peak signal
% into its component parts. The objective is to determine whether your
% signal can be represented as the sum of fundamental underlying peaks
% shapes. Accepts signals of any length, including those with non-integer
% and non-uniform x-values. Fits any number of peaks of any of 45 curve
% shapes.
%
% Version 8.2: August 2016, adds shape 48 equal-width exponential pulse.
%
% For more details, see
% http://terpconnect.umd.edu/~toh/spectrum/CurveFittingC.html and
% http://terpconnect.umd.edu/~toh/spectrum/InteractivePeakFitter.htm
% 
global AA xxx PEAKHEIGHTS FIXEDPARAMETERS AUTOZERO delta BIPOLAR CLIPHEIGHT
% format short g
format compact
warning off all
NumArgOut=nargout;
datasize=size(signal);
if datasize(1)<datasize(2),signal=signal';end
datasize=size(signal);
if datasize(2)==1, %  'signal' is a vector; Must be peakfit(Y-vector)
    X=1:length(signal); % Create an independent variable vector
    Y=signal;
else
    % 'signal' is a matrix. Must be peakfit(DataMatrix)
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
        extra=1;
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
        extra=1;
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
        extra=1;
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
        extra=1;
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
        extra=ones(1,NumPeaks);
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
        % initialstart=start % testing
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
shapesvector=peakshape;
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
firststart=start; % <<<<<<<<<<<
if FIXEDPARAMETERS==0, FIXEDPARAMETERS=length(xx)/10;end
if peakshape==16;FIXEDPOSITIONS=fixedparameters;end
if peakshape==17;FIXEDPOSITIONS=fixedparameters;end
if AUTOZERO>3,AUTOZERO=3;disp('AUTPOZERO must be between 0 and 3');end
if AUTOZERO<0,AUTOZERO=0;disp('AUTPOZERO must be between 0 and 3');end
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
        ShapeString='Lognormal';
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
    case 34
        ShapeString='Fixed-width Voigt';
    case 35
        ShapeString='Fixed-width G/L blend';
    case 36
        ShapeString='Fixed-width ExpGaussian';
    case 37
        ShapeString='Fixed-width Pearson';
    case 38
        ShapeString='ExpLorentzian (var. time constant)'; 
    case 40
        ShapeString='Sine wave';
    case 41
        ShapeString='Rectangular pulse';
    case 42
        ShapeString='Flattened Gaussian';  
    case 43
        ShapeString='3-parameter Gompertz.';  
    case 44
        ShapeString='1-exp(-k*t)';  
    case 45
        ShapeString='4-parameter logistic'; 
    case 46
        ShapeString='Quadratic Baseline'; 
    case 47
        ShapeString='Blackbody'; 
    case 48
        ShapeString='Equal-width exp. pulse'; 
    otherwise
end % switch peakshape
  
% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.0000001,'Display','off','MaxFunEvals',1000 );
LowestError=1000; % or any big number greater than largest error expected
FitParameters=zeros(1,NumPeaks.*3); 
BestStart=zeros(1,NumPeaks.*2); 
height=zeros(1,NumPeaks); 
bestmodel=zeros(size(yy));
TrialParameters=zeros(1,NumPeaks.*3);
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
            % fixedstart=fixedstart
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
            TrialParameters=fminsearch(@(lambda)(fitmultiple(lambda,xx,yy,NumPeaks,shapesvector,extra)),newstart,options);
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
            % newstart=newstart % testing
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
%             case30newstart=newstart % uncomment for testing
%             sizestartcase30=size(start) % uncomment for testing

            TrialParameters=fminsearch(@(lambda)(fitvoigtv(lambda,xx,yy)),start);
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
            % case31newstart=newstart % uncomment for testing
            % sizestartcase31=size(start) % uncomment for testing
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
            case32newstart=newstart % uncomment for testing
            sizestartcase32=size(start) % uncomment for testing
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
            % newstart=newstart % uncomment for testing
            TrialParameters=fminsearch(@(lambda)(fitGLv(lambda,xx,yy)),newstart);
        case 34
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end  
            % fixedstart=fixedstart % uncomment for testing
            TrialParameters=fminsearch(@(lambda)(fitFWVoigt(lambda,xx,yy,extra)),fixedstart,options);
        case 35
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end            
            TrialParameters=fminsearch(@(lambda)(fitFWGL(lambda,xx,yy,extra)),fixedstart,options);
        case 36
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end            
            TrialParameters=fminsearch(@(lambda)(fitFWExpGaussian(lambda,xx,yy,extra)),fixedstart,options);
        case 37
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end            
            TrialParameters=fminsearch(@(lambda)(fitFWPearson(lambda,xx,yy,extra)),fixedstart,options);
        case 38
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
            TrialParameters=fminsearch(@(lambda)(fitexplorentzianv(lambda,zxx,zyy)),newstart);
        case 40
             TrialParameters=fminsearch(@(lambda)(fitsine(lambda,xx,yy)),newstart,options); 
        case 41
            TrialParameters=fminsearch(@(lambda)(fitrectangle(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 42
            TrialParameters=fminsearch(@(lambda)(fitngaussian(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 43
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
            TrialParameters=fminsearch(@(lambda)(fitGompertz(lambda,xx,yy)),newstart);

        case 44
            TrialParameters=fminsearch(@(lambda)(fitOneMinusExp(lambda,xx,yy)),newstart,options);
         case 45
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
            % newstart=newstart % uncomment for testing
            TrialParameters=fminsearch(@(lambda)(fitFourPL(lambda,xx,yy)),newstart);
        case 46
            TrialParameters=fminsearch(@(lambda)(fitquadslope(lambda,xx,yy)),newstart,options);
        case 47
            bbstart=3000;
            TrialParameters=fminsearch(@(lambda)(fitblackbody(lambda,xx,yy)),bbstart,options);
        case 48
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewexpulse(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        otherwise
    end % switch peakshape

    % Check variables
    % sizeNewstart=size(newstart) % uncomment for testing
  
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
            A(m,:)=gaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
        case 12
            A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
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
        case 34
             width(m)=abs(FIXEDPARAMETERS(m));
%                 gD(m)=width(m);
%                 gL(m)=extra.*gD(m);
%                 width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2))
            A(m,:)=voigt(xx,TrialParameters(m), width(m),extra);
        case 35
            A(m,:)=GL(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);    
        case 36
            A(m,:)=expgaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);    
        case 37
            A(m,:)=pearson(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);    
        case 38
            A(m,:)=explorentzian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),-TrialParameters(3*m));        
        case 40
            A(m,:)=sine(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 41
            A(m,:)=rectangle(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 42
            A(m,:)=ngaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 43
            A(m,:)=Gompertz(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
        case 44
            A(m,:)=OneMinusExp(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 45
            A(m,:)=FourPL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
        case 46
            A(m,:)=quadslope(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 47
            A(m,:)=blackbody(xx,TrialParameters(m));
        case 48
            A(m,:)=exppulse(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
    end % switch
    xxrange=max(xx)-min(xx);
    for parameter=1:2:2*NumPeaks,
        newstart(parameter)=newstart(parameter)+(xxrange.*(delta*randn)./(NumPeaks+1));
        newstart(parameter+1)=newstart(parameter+1)*(1+delta*(rand-.5)/100);
    end
end % for NumPeaks
% newstart=newstart; % <<<<<<<<<< % error check
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
%          sizePeakHeights=size(PEAKHEIGHTS) %  % uncomment for testing
%          SizeA=size(A)  % uncomment for testing
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
%  ErrorVector(k)=MeanFitError; %  % uncomment for testing
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
        AA(m,:)=gaussian(xxx,FitParameters(m),FIXEDPARAMETERS(m));
    case 12
        AA(m,:)=lorentzian(xxx,FitParameters(m),FIXEDPARAMETERS(m));
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
    case 34
                  width(m)=abs(FIXEDPARAMETERS(m));
%                 gD(m)=width(m);
%                 gL(m)=extra.*gD(m);
%                 width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 +
%                 gD(m).^2));
        AA(m,:)=voigt(xxx,FitParameters(m),width(m),extra);
    case 35
        AA(m,:)=GL(xxx,FitParameters(m),FIXEDPARAMETERS(m),extra);
    case 36
        AA(m,:)=expgaussian(xxx,FitParameters(m),FIXEDPARAMETERS(m),extra);    
    case 37
        AA(m,:)=pearson(xxx,FitParameters(m),FIXEDPARAMETERS(m),extra);
    case 38
        AA(m,:)=explorentzian(xxx,FitParameters(3*m-2),FitParameters(3*m-1),-FitParameters(3*m)*length(xxx)./length(xx));        
    case 40
        AA(m,:)=sine(xx,FitParameters(2*m-1),FitParameters(2*m));
    case 41
        AA(m,:)=rectangle(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 42
        AA(m,:)=ngaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 43
        AA(m,:)=Gompertz(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m));  
    case 44
        AA(m,:)=OneMinusExp(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 45
        AA(m,:)=FourPL(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m));    
    case 46
        AA(m,:)=quadslope(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 47
        AA(m,:)=blackbody(xxx,FitParameters(m));
    case 48
        AA(m,:)=exppulse(xxx,FitParameters(m),FitParameters(NumPeaks+1));
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
residual=yy-bestmodel;

if plots,
    % Mark starting peak positions with vertical dashed magenta lines
    if peakshape(1)==16||peakshape(1)==17
    else
        if peakshape(1)==29, % Segmented linear
%            subplot(2,1,1);plot([PEAKHEIGHTS' PEAKHEIGHTS'],[0 max(yy)],'m--')
        else
            for marker=1:NumPeaks,
                markx=BestStart((2*marker)-1);
%                subplot(2,1,1);plot([markx+xoffset markx+xoffset],[0 max(yy)],'m--')
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
            title(['ipf 12.2   No baseline correction'])
        case 1
            title(['ipf 12.2   Linear baseline subtraction'])
        case 2
            title(['ipf 12.2   Quadratic subtraction baseline'])
        case 3
            title(['ipf 12.2   Flat baseline correction'])
    end
 
    switch peakshape(1)
        case {4,20,34,37,42}
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      Shape Constant = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000) ] )
        case {5,8,18,36}
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      Time Constant = ' num2str(extra)   '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000)  ] )
        case 13
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      % Gaussian = ' num2str(extra)   '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000)  ] )
        case {14,15,22,35}
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      extra = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000) ] )
        case 28
            xlabel(['Shape = ' ShapeString '      Order = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '%  R2 = ' num2str(round(1000*LowestError)/1000) ] )
        case 43
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '        Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000) ] )
         otherwise
            if peakshape(1)==29, % Segmented linear
                xlabel(['Breakpoints = ' num2str(NumPeaks) '     Shape = ' ShapeString  '     Error = ' num2str(round(1000*LowestError)/1000) '%  R2 = ' num2str(round(100000*Rsquared)/100000) ] )
            else
                xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH)  '     Error = ' num2str(round(1000*LowestError)/1000) '%  R2 = ' num2str(round(100000*Rsquared)/100000) ] )
            end % if peakshape(1)==29
    end % switch peakshape(1)

    % Bottom half of the figure shows the residuals and displays RMS error
    % between original signal and model
    % residual=yy-bestmodel;
    subplot(2,1,2);plot(xx+xoffset,residual,'m.')
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
%  FitParameters=FitParameters
switch peakshape(1),
    case {6,7,8,48}, % equal-width peak models only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
            end
        end
    case {11,12,34,35,36,37}, % Fixed-width shapes only
        for m=1:NumPeaks,
            width(m)=abs(FitParameters(m));
            if peakshape==34,
                gD(m)=width(m);
                gL(m)=extra.*gD(m);
                width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2));
            end
            if m==1,
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
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
    case {30,31,32,33,38,43} % Special case of shapes with 3 iterated variables
        for m=1:NumPeaks,
            width(m)=abs(FitParameters(3*m-1));
            if peakshape==30,
                gD(m)=width(m);
                gL(m)=FitParameters(3*m).*gD(m);
                width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2));
            end
            if m==1,
                FitResults=[round(m) FitParameters(3*m-2) height(m) width(m) area(m) FitParameters(3*m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(3*m-2) height(m) width(m) area(m)] FitParameters(3*m)];
            end
        end
    case 47 % Shapes with 1 iterated variable
        FitParameters=FitParameters;
    otherwise % Normal shapes with 2 iterated variables
        for m=1:NumPeaks,
            width(m)=abs(FitParameters(2*m));
            if peakshape==20,
                gD=width(m);
                gL=extra.*gD;
                width(m) = 2.*(0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2));
            end
            if m==1,
                FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)]];
            end % if m==1
            
        end % for m=1:NumPeaks,
end % switch peakshape(1)
  
% Rearrange fit results for Gompertz to Bo, Kh, and L
if peakshape(1)==43;
    for m=1:NumPeaks,
        FitResults(m,2)=FitResults(m,2).*FitResults(m,3);
        FitResults(m,4)=FitResults(m,3).*FitResults(m,4);
        FitResults(m,3)=1;
    end
end

%Sort FitResults
FitResults=sortrows(FitResults,2);

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
        case {9,19,10,23,40}  % Pulse and sigmoid shapes only
            text(startx,starty+dyy/2,['Peak #          tau1           Height           tau2             Area'] );
        case 28, % Polynomial
            text(startx,starty+dyy/2,['Polynomial coefficients'] );
        case 29 % Segmented linear
             text(startx,starty+dyy/2,['x-axis breakpoints'] );
        case {30,31,32,33,38} % Special case of shapes with 3 iterated variables
            text(startx,starty+dyy/2,['Peak #          Position        Height         Width             Area       Shape factor'] );            
        case 43 % 3 parameter Gompertz
            text(startx,starty+dyy/2,['Peak #           Bo             Height            Kh                Area                 L'] );            
        otherwise
            text(startx,starty+dyy/2,['Peak #          Position         Height         Width             Area   '] );
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
    if peakshape(1)==30 || peakshape(1)==31 || peakshape(1)==32 || peakshape(1)==33 || peakshape(1)==38 || peakshape(1)==43,
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
        [FitResults,BootFitError]=fitpeaks(bx,by,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,FIXEDPARAMETERS,shapesvector);
        for peak=1:NumPeaks,
            switch peakshape(1)
                case {30,31,32,33,38,43}
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
            % Label columns for bootstrap results
            switch peakshape(1)
                case {9,19,10,23,40}  % Pulse and sigmoid shapes only
                    disp(['Peak #',num2str(peak) '         tau1           Height         tau2           Area'] );
                case {30,31,32,33,38} % Special case of shapes with 3 iterated variables
                    disp(['Peak #',num2str(peak) '         Position        Height         Width             Area       Shape factor'] );
                case 43 % 3 parameter Gompertz
                    disp(['Peak #',num2str(peak) '         Bo             Height         Kh                       L'] );
                otherwise
                    disp(['Peak #',num2str(peak) '         Position         Height         Width             Area'] );
            end
            
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
        format short g
        if plots,
            disp(['Mean:       ', num2str(BootstrapMean)])
            disp(['STD:         ', num2str(BootstrapSTD)])
            disp(['STD (IQR):   ', num2str(BootstrapIQR)])
            disp(['% RSD:       ', num2str(PercentRSD)])
            disp(['% RSD (IQR): ', num2str(PercentIQR)])
        end % if plots
        BootResults(peak,:)=[BootstrapMean BootstrapSTD PercentRSD BootstrapIQR PercentIQR];
    end % peak=1:NumPeaks,
end % if NumArgOut==8,
if AUTOZERO==3;
else
    baseline=bkgcoef;
end
% ----------------------------------------------------------------------
function [FitResults,LowestError]=fitpeaks(xx,yy,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,fixedparameters,shapesvector)
% Based on peakfit Version 3: June, 2012. 
global PEAKHEIGHTS FIXEDPARAMETERS BIPOLAR MINWIDTH coeff
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
options = optimset('TolX',.000001,'Display','off','MaxFunEvals',1000 );
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
            TrialParameters=fminsearch(@(lambda)(fitupsigmoid(lambda,xx,yy)),newstart,options);
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
            TrialParameters=fminsearch(@(lambda)(fitmultiple(lambda,xx,yy,NumPeaks,shapesvector,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 23
            TrialParameters=fminsearch(@(lambda)(fitdownsigmoid(lambda,xx,yy)),newstart,optionst);
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
        case 34
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(fitFWVoigt(lambda,xx,yy,extra)),fixedstart,options);
        case 35
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(fitFWGL(lambda,xx,yy,extra)),fixedstart,options);
        case 36
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(fitFWExpGaussian(lambda,xx,yy,extra)),fixedstart,options);
        case 37
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(fitFWPearson(lambda,xx,yy,extra)),fixedstart,options);
        case 38
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexplorentzianv(lambda,zxx,zyy)),newstart);         
        case 40
             TrialParameters=fminsearch(@(lambda)(fitsine(lambda,xx,yy)),newstart,options);
        case 41
            TrialParameters=fminsearch(@(lambda)(fitrectangle(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 42
            TrialParameters=fminsearch(@(lambda)(fitngaussian(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 43
            TrialParameters=fminsearch(@(lambda)(fitGompertz(lambda,xx,yy)),newstart);
        case 44
            TrialParameters=fminsearch(@(lambda)(fitOneMinusExp(lambda,xx,yy)),newstart);           
        case 45
            TrialParameters=fminsearch(@(lambda)(fitFourPL(lambda,xx,yy)),newstart);  
        case 46
            TrialParameters=fminsearch(@(lambda)(fitquadslope(lambda,xx,yy)),newstart,options);
        case 47
            bbstart=3000;
            TrialParameters=fminsearch(@(lambda)(fitblackbody(lambda,xx,yy)),bbstart,options);
        case 48
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewexpulse(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
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
                A(m,:)=gaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
            case 12
                A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
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
            case 34
                width(m)=abs(FIXEDPARAMETERS(m));

%                 gD(m)=width(m);
%                 gL(m)=extra.*gD(m);
%                 width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2))

                A(m,:)=voigt(xx,TrialParameters(m), width(m),extra);
            case 35
                A(m,:)=GL(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);
            case 36
                A(m,:)=expgaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);
            case 37
                A(m,:)=pearson(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);
            case 38
                A(m,:)=explorentzian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));                    
            case 40
                A(m,:)=sine(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 41
                A(m,:)=rectangle(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 42
                A(m,:)=ngaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 43
                A(m,:)=Gompertz(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));
            case 44
                A(m,:)=OneMinusExp(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 45
                A(m,:)=FourPL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
            case 46
                A(m,:)=quadslope(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 47
                A(m,:)=blackbody(xx,TrialParameters(m));
            case 48
                A(m,:)=exppulse(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
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
    case {11,12,35,36,37}, % Fixed-width shapes only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
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
    case {30,31,32,33,38,43} % Special case of shapes with 3 iterated variables
        for m=1:NumPeaks,
            width(m)=abs(FitParameters(3*m-1));
            if peakshape==30,
                gD(m)=width(m);
                gL(m)=FitParameters(3*m).*gD(m);
                width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2));
            end
            if m==1,
                FitResults=[round(m) FitParameters(3*m-2) height(m) width(m) area(m) FitParameters(3*m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(3*m-2) height(m) width(m) area(m) FitParameters(3*m)]];
            end
        end
    otherwise % Normal shapes with 2 iterated variables
        for m=1:NumPeaks,
            width(m)=abs(FitParameters(2*m));
            if peakshape==20,
                gD=width(m);
                gL=extra.*gD;
                width(m) = 2.*(0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2));
            end
            if m==1,
                FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)]];
            end % if m==1
        end % for m=1:NumPeaks,
end % switch peakshape(1)
if peakshape==34,
        DW=2*(0.5346*a*1.2772 + sqrt(0.2166*a*1.2772.^2 + 1.2772.^2))
end
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
% Fitting function for multiple Gaussian peaks.
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
function err = fitngaussian(lambda,t,y,shapeconstant)
%   Fitting functions for multiple flattened Gaussian peaks.
% T. C. O'Haver (toh@umd.edu),   Version 1.3, October 23, 2006.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = ngaussian(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
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
% Fitting function for multiple Gaussian peaks with equal peak widths.
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
%	Fitting function for multiple fixed width Gaussians
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),FIXEDPARAMETERS(j))';
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
%	Fitting function for multiple fixed-position Gaussians
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
%	Fitting function for multiple fixed-position Lorentzians
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
%	Fitting function for multiple fixed width Lorentzians
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),FIXEDPARAMETERS(j))';
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
% Fitting function for multiple  Lorentzian band signal with equal peak widths.
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
g = exp(-((x-pos)./(0.60056120439323.*wid)).^2);
% ----------------------------------------------------------------------
function g = ngaussian(x,pos,wid,n)
%  ngaussian(x,pos,wid) = flattened Gaussian centered on x=pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
% Shape is Gaussian when n=1. Becomes more rectangular as n increases.
%  T. C. O'Haver, 1988, revised 2014
% Example: ngaussian([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
if n>0,
    g = 1-(10.^-(n.*gaussian(x,pos,wid)));
    g=g./max(g);
else
    g = gaussian(x,pos,wid);
end
% ----------------------------------------------------------------------
function err = fitlorentzian(lambda,t,y)
%	Fitting function for multiple lorentzians, lambda(1)=position, lambda(2)=width
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
% [lambda PEAKHEIGHTS err]
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
%	Fitting function for multiple logistic peaks, lambda(1)=position, lambda(2)=width
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
%	Fitting function for multiple triangular, lambda(1)=position, lambda(2)=width
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
%Triangle function.  pos=position; wid=half-width (both scalar)
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
function err = fitrectangle(lambda,t,y)
%	Fitting function for multiple rectangle, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fitrectangle assumes a rectangle function 
%  T. C. O'Haver, May 2016
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = rectangle(t,lambda(2*j-1),lambda(2*j))';
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
function g = rectangle(x,pos,wid)
%rectangle function.  pos=position; wid=half-width (both scalar)
%rectangle(x,pos,wid), where x may be scalar or vector,
%pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 2016
% Example
% x=[0:.1:10];plot(x,rectangle(x,5.5,2.3),'.')
g=zeros(size(x));
hw=wid./2;
for i=1:length(x),  
if x(i)<pos-hw,g(i)=0;end
if x(i)>pos-hw,g(i)=1;end
if x(i)>pos+hw,g(i)=0;end
end
% ----------------------------------------------------------------------
function err = fitpearson(lambda,t,y,shapeconstant)
%   Fitting functions for multiple Pearson 7 bands.
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
% Fitting functions for multiple pearson functions with independently variable
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
function err = fitFWPearson(lambda,t,y,shapeconstant)
%	Fitting function for multiple fixed width Pearson7
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = pearson(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
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
% g = pearson(x,pos,wid,m) where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% m=some number
%  T. C. O'Haver, 1990  
g=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m;
% ----------------------------------------------------------------------
function err = fitexpgaussian(lambda,t,y,timeconstant)
%   Fitting functions for multiple exponentially-broadened Gaussian bands signal.
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
%   Fitting functions for multiple exponentially-broadened lorentzian band signal.
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
function err = fitexplorentzianv(lambda,t,y)
% Fitting functions for multiple exponentially-broadened Gaussians with
% independently variable time constants
% T. C. O'Haver (toh@umd.edu), 2015.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = explorentzian(t,lambda(3*j-2),lambda(3*j-1),-lambda(3*j))';
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
% Fitting function for multiple exponentially-broadened Gaussian bands with equal peak widths.
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
% Fitting functions for multiple exponentially-broadened Gaussians with
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
function err = fitFWExpGaussian(lambda,t,y,shapeconstant)
%	Fitting function for multiple fixed width exponentially-broadened gaussian
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = expgaussian(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
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
g = exp(-((x-pos)./(0.60056120439323.*wid)) .^2);
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
function err = fitewexpulse(lambda,t,y)
% Fitting function for multiple Gaussian peaks with equal peak widths.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = exppulse(t,lambda(j),lambda(numpeaks+1))';
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
% Fitting function for iterative fit to the sum of multiple
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
% Fitting function for iterative fit to the sum of multiple
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
%   Fitting functions for multiple Gaussian/Lorentzian blend peaks
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
function err = fitFWGL(lambda,t,y,shapeconstant)
%	Fitting function for a multiple fixed width Gaussian/Lorentzian blend
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = GL(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
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
% Fitting functions for multiple Gaussian/Lorentzian blend functions with
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
% Fitting functions for multiple Voigt profile function
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
function err = fitFWVoigt(lambda,t,y,shapeconstant)
%	Fitting function for multiple fixed width Voigt
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
% numpeaksfitFWVoigt=numpeaks
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = voigt(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
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
% Fitting functions for multiple Voigt profile function with independently variable
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
% sizeabs=size(abs(xx-pos))
% sizegV=size(gV)
y = abs(xx-pos)./gV;
g = 1/(2*gV*(1.065 + 0.447*x + 0.058*x^2))*((1-x)*exp(-0.693.*y.^2) + (x./(1+y.^2)) + 0.016*(1-x)*x*(exp(-0.0841.*y.^2.25)-1./(1 + 0.021.*y.^2.25)));
g=g./max(g);
% ----------------------------------------------------------------------
function err = fitBiGaussian(lambda,t,y,shapeconstant)
%   Fitting functions for multiple BiGaussian.
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
% BiGaussian (different width on leading edge and trailing edge).
% pos=position; wid=width 
% m = ratio of widths of right-hand to left-hand halves.
% If m=1, it becomes identical to a Gaussian.
% Verison 2, T. C. O'Haver, 2012
% Example: Plots Gaussian and BiGaussian with m=3, both with halfwidth=20.
% x=[1:100];
% g=gaussian(x,50,20);
% bg=bigaussian(x,50,20,3);
% plot(x,g,x,bg)
%
lx=length(x);
hx=val2ind(x,pos);
hwid=2.*wid;
g(1:hx)=gaussian(x(1:hx),pos,hwid/(m+1));
g(hx+1:lx)=gaussian(x(hx+1:lx),pos,m.*hwid/(m+1));
% ----------------------------------------------------------------------
function err = fitBWF(lambda,t,y,shapeconstant)
%   Fitting function for multiple Breit-Wigner-Fano.
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
% Fitting function for iterative fit to the sum of multiple
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
% Fitting function for iterative fit to the sum of multiple
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
% Fitting function for iterative fit to the sum of multiple
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
function g=sine(t,f,phase) 
% Sine wave (alpha test)
g=sin(2*pi*f*(t+phase));
% ----------------------------------------------------------------------
function err = fitd1gauss(lambda,t,y)
%   Fitting functions for multiple first derivative of a Gaussian
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
yi=interp1(segs,yy,x);
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
function err = fitquadslope(lambda,x,y)
% Fitting function for iterative fit to linear function
global PEAKHEIGHTS LOGPLOT
A = zeros(length(x),round(length(lambda)/2));
for j = 1:length(lambda)/2,
     A(:,j)=quadslope(x,lambda(1),lambda(2));
end
PEAKHEIGHTS=A\y';
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function y=quadslope(x,boa,coa) % normalized quadratic
 y=(x.^2+(boa).*x+coa);
y=y./max(y);
% ----------------------------------------------------------------------
function err = fitGompertz(lambda,t,y)
% Fitting functions for multiple Gompertz function
% T. C. O'Haver (toh@umd.edu), 2016.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
% Sizelambda=size(lambda)
for j = 1:length(lambda)/3,
    A(:,j) = Gompertz(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
% PEAKHEIGHTS=1;
% SizeA=size(A)
% sizwPEAKHEIGHTS=size(PEAKHEIGHTS)
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function y=Gompertz(t,Bo,Kh,L)
% A Gompertz curve or Gompertz function, named after Benjamin Gompertz, is
% a sigmoid function. It is a type of mathematical model for a time series,
% where growth is slowest at the start and end of a time period. The
% right-hand or future value asymptote of the function is approached much
% more gradually by the curve than the left-hand or lower valued asymptote,
% in contrast to the simple logistic function in which both asymptotes are
% approached by the curve symmetrically. It is a special case of the
% generalized logistic function.
% 
% Example:
% x=1:.1:10;y=gompertz(x,6,3,4);plot(x,y)
y=Bo*exp(-exp((Kh*exp(1)/Bo)*(L-t) +1));
% ----------------------------------------------------------------------
function err = fitFourPL(lambda,t,y)
% Fitting functions for multiple Gompertz function
% T. C. O'Haver (toh@umd.edu), 2016.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = FourPL(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
PEAKHEIGHTS=1;
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function y=FourPL(x,miny,slope,ip)
% Normalized four parameters logistic regression  
% https://psg.hitachi-solutions.com/masterplex/blog/the-4-parameter-logisti
% c-4pl-nonlinear-regression-model
% miny = minimum asymptote. In an ELISA assay where you have a standard
% curve, this can be thought of as the response value at 0 standard
% concentration. slope = Hill slope. The Hill Slope or slope factor refers
% to the steepness of the curve. It could either be positive or negative.
% As the absolute value of the Hill slope increases, so does the steepness
% of the curve. ip = inflection point: The inflection point is defined as the
% point on the curve where the curvature changes direction or signs. This
% can be better explained if you can imagine the concavity of a sigmoidal
% curve. The inflection point is where the curve changes from being concave
% upwards to concave downwards. maxy = maximum asymptote. In an ELISA assay
% where you have a standard curve, this can be thought of as the response
% value for infinite standard concentration.
% Example:
% x=0:20;
% miny=0;slope=5;ip=10;d=0;maxy=10;
% y=FourPL(x,miny,slope,ip,maxy);plot(x,y)
%
y = 1+(miny-1)./(1+(x./ip).^slope);
% ----------------------------------------------------------------------
function err = fitOneMinusExp(lambda,t,y)
%   Fitting functions for multiple first derivative of a Gaussian
%  T. C. O'Haver, 2014
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = OneMinusExp(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = OneMinusExp(x,pos,wid)
% OneMinusExp(x,pos,wid) = 1-exp(-wid.*(x-pos));
% Example:  x=0:10;y=1-exp(-(.5.*x));plot(x,y)
g = 1-exp(-wid.*(x-pos));
% ----------------------------------------------------------------------
function err = fitblackbody(lambda,wavelength,y)
%  Fitting function for a blackbody spectrum.
%  T. C. O'Haver, May 2008
global PEAKHEIGHTS AUTOZERO BIPOLAR
% sizelambda=size(lambda)
radiance = blackbody(wavelength,lambda(1));
% if AUTOZERO==3,A=[ones(size(y))' A];end
% if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
% sizePEAKHEIGHTS=size(PEAKHEIGHTS)
% sizey=size(y)
PEAKHEIGHTS = radiance/y;
z = radiance*PEAKHEIGHTS;
err = norm(z-y);
% ----------------------------------------------------------------------
function radiance=blackbody(wavelength,temperature)
radiance = 1.19111E+16*wavelength.^(-5)./(exp(14380000./(wavelength*temperature))-1);
% ----------------------------------------------------------------------
function b=iqr(a)
% b = IQR(a) returns the interquartile range divided by 1.34896, of the
% values in a. If a is a vector, b is the difference between the 75th and
% 25th percentiles of a, divided by 1.34896. If b is a matrix, b is a row
% vector containing the interquartile range of each column of a, divided by
% 1.34896. The factor 1.34896 makes this quantity equal on average to the
% standard deviation (SD) of normal distributions, and thue the IRQ is a
% better estimate of the standard deviation without ouliers.
%  T. C. O'Haver, 2017
mina=min(a);
sizea=size(a);
NumCols=sizea(2);
for n=1:NumCols,b(:,n)=a(:,n)-mina(n);end
Sorteda=sort(b);
lx=length(Sorteda);
SecondQuartile=round(lx/4);
FourthQuartile=3*round(lx/4);
b=abs(Sorteda(FourthQuartile,:)-Sorteda(SecondQuartile,:));
b=b./1.34896;
% ----------------------------------------------------------------------
function err = fitmultiple(lambda,xx,y,numpeaks,shapesvector,extra)
% Fitting function for a multiple-shape band signal.
% The sequence of peak shapes are defined by the vector "shape".
% The vector "extra" determines the shape of variable-shape peaks.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT coeff FIXEDPARAMETERS
% FIXEDPARAMETERS=FIXEDPARAMETERS % testing
% PEAKHEIGHTS=PEAKHEIGHTS % testing
% lambda=lambda
% numpeaks=round(length(lambda)/2))
A=zeros(numpeaks,length(xx));
for j = 1:numpeaks,
    if shapesvector(j)==28,
        coeff=polyfit(xx,y,extra(j));
        A(j,:) = polyval(coeff,xx);
    else
        % sizeA=size(A)
        switch shapesvector(j)
            case 1
                A(j,:)=gaussian(xx,lambda(2*j-1),lambda(2*j));
            case 2
                A(j,:)=lorentzian(xx,lambda(2*j-1),lambda(2*j));
            case 3
                A(j,:)=logistic(xx,lambda(2*j-1),lambda(2*j));
            case 4
                A(j,:)=pearson(xx,lambda(2*j-1),lambda(2*j),extra(j));
            case 5
                A(j,:)=expgaussian(xx,lambda(2*j-1),lambda(2*j),-extra(j))';
            case 6
                A(j,:)=gaussian(xx,lambda(j),lambda(NumPeaks+1));
            case 7
                A(j,:)=lorentzian(xx,lambda(j),lambda(NumPeaks+1));
            case 8
                A(j,:)=expgaussian(xx,lambda(j),lambda(NumPeaks+1),-extra(j))';
            case 9
                A(j,:)=exppulse(xx,lambda(2*j-1),lambda(2*j));
            case 10
                A(j,:)=upsigmoid(xx,lambda(2*j-1),lambda(2*j));
            case 11
                A(j,:)=gaussian(xx,lambda(j),FIXEDPARAMETERS(j));
            case 12
                A(j,:)=lorentzian(xx,lambda(j),FIXEDPARAMETERS(j));
            case 13
                A(j,:)=GL(xx,lambda(2*j-1),lambda(2*j),extra(j));
            case 14
                A(j,:)=BiGaussian(xx,lambda(2*j-1),lambda(2*j),extra(j));
            case 15
                A(j,:)=BWF(xx,lambda(2*j-1),lambda(2*j),extra(j));
            case 16
                A(j,:)=gaussian(xx,FIXEDPOSITIONS(j),lambda(j));
            case 17
                A(j,:)=lorentzian(xx,FIXEDPOSITIONS(j),lambda(j));
            case 18
                A(j,:)=explorentzian(xx,lambda(2*j-1),lambda(2*j),-extra(j))';
            case 19
                A(j,:)=alphafunction(xx,lambda(2*j-1),lambda(2*j));
            case 20
                A(j,:)=voigt(xx,lambda(2*j-1),lambda(2*j),extra(j));
            case 21
                A(j,:)=triangular(xx,lambda(2*j-1),lambda(2*j));
            case 22
                A(j,:)=peakfunction(shapesvector(j),xx,lambda(2*j-1),lambda(2*j),extra(j));
            case 23
                A(j,:)=downsigmoid(xx,lambda(2*j-1),lambda(2*j));      
            case 24
                A(j,:)=nbinpdf(xx,lambda(2*j-1),lambda(2*j));
            case 25
                A(j,:)=lognormal(xx,lambda(2*j-1),lambda(2*j));
            case 26
                A(j,:)=linslope(xx,lambda(2*j-1),lambda(2*j));
            case 27
                A(j,:)=d1gauss(xx,lambda(2*j-1),lambda(2*j));       
            case 28
                A(j,:)=polynomial(xx,lambda(2*j-1),lambda(2*j));       
            case 29
                A(j,:)=segmented(xx,yy,PEAKHEIGHTS);
            case 30
                A(j,:)=voigt(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));        
            case 31
                A(j,:)=expgaussian(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));        
            case 32
                A(j,:)=pearson(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));        
            case 33
                A(j,:)=GL(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));
            case 34
                width(j)=abs(FIXEDPARAMETERS(j));
                 % [j lambda(j) width(j) extra(j)]
%                 gD(j)=width(j);
%                 gL(j)=extra.*gD(j);
%                 width(j) = 2.*(0.5346*gL(j) + sqrt(0.2166*gL(j).^2 +
%                 gD(j).^2))
                A(j,:)=voigt(xx,lambda(j),width(j),extra(j));
                % figure(5);plot(A(j,:));figure(1)
            case 35
                A(j,:)=GL(xx,lambda(j),FIXEDPARAMETERS(j),extra(j));
            case 36
                A(j,:)=expgaussian(xx,lambda(j),FIXEDPARAMETERS(j),extra(j));
            case 37
                A(j,:)=pearson(xx,lambda(j),FIXEDPARAMETERS(j),extra(j));
            case 38
                A(j,:)=explorentzian(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));                    
            case 40
                A(j,:)=sine(xx,lambda(2*j-1),lambda(2*j));
            case 41
                A(j,:)=rectangle(xx,lambda(2*    j-1),lambda(2*j));
            case 42
                A(j,:)=ngaussian(xx,lambda(2*j-1),lambda(2*j),extra(j));
            case 43
                A(j,:)=Gompertz(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));
            case 44
                A(j,:)=OneMinusExp(xx,lambda(2*j-1),lambda(2*j));
            case 45
                A(j,:)=FourPL(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));        
            case 46
                A(j,:)=quadslope(xx,lambda(2*j-1),lambda(2*j));
            case 47
                A(j,:)=blackbody(xx,lambda(j));
            case 48
                A(j,:)=exppulse(xx,lambda(j),lambda(NumPeaks+1));
        end % switch
    end % if shapesvector
end % for j=1:numpeaks,
if AUTOZERO==3,A=[ones(size(y))' A];end
% sizeA=size(A)
% max(A')
% sizeyp=size(y')
if BIPOLAR,PEAKHEIGHTS=A'\y';else PEAKHEIGHTS=abs(A'\y');end
% PEAKHEIGHTS=PEAKHEIGHTS
z = A'*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function p=peakfunction(shape,x,pos,wid,extra,coeff)
global FIXEDPARAMETERS
% function that generates any of 20 peak types specified by number. 'shape'
% specifies the shape type of each peak in the signal: "peakshape" = 1-20.
% 1=Gaussian 2=Lorentzian, 3=logistic, 4=Pearson, 5=exponentionally
% broadened Gaussian; 9=exponential pulse, 10=up sigmoid,
% 13=Gaussian/Lorentzian blend; 14=BiGaussian, 15=Breit-Wigner-Fano (BWF) ,
% 18=exponentionally broadened Lorentzian; 19=alpha function; 20=Voigt
% profile; 21=triangular; 23=down sigmoid; 25=lognormal. "extra" is required
% for variable-shape peaks only.
switch shape,
    case 1
        p=gaussian(x,pos,wid);
    case 2
        p=lorentzian(x,pos,wid);
    case 3
        p=logistic(x,pos,wid);
    case 4
        p=pearson(x,pos,wid,extra);
    case 5
        p=expgaussian(x,pos,wid,extra);
    case 6
        p=gaussian(x,pos,wid);
    case 7
        p=lorentzian(x,pos,wid);
    case 8
        p=expgaussian(x,pos,wid,extra)';
    case 9
        p=exppulse(x,pos,wid);
    case 10
        p=upsigmoid(x,pos,wid);
    case 11
        p=gaussian(x,pos,wid);
    case 12
        p=lorentzian(x,pos,wid);
    case 13
        p=GL(x,pos,wid,extra);
    case 14
        p=BiGaussian(x,pos,wid,extra);
    case 15
        p=BWF(x,pos,wid,extra);
    case 16
        p=gaussian(x,pos,wid);
    case 17
        p=lorentzian(x,pos,wid);
    case 18
        p=explorentzian(x,pos,wid,extra)';
    case 19
        p=alphafunction(x,pos,wid);
    case 20
        p=voigt(x,pos,wid,extra);
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
    case 30
        p=voigt(x,pos,pos,wid,extra);
    case 31
        p=expgaussian(x,pos,wid,-extra);
    case 32
        p=pearson(x,pos,wid,extra);
    case 33
        p=GL(x,pos,wid,extra);
    case 34
        % width(extra)=abs(FIXEDPARAMETERS(extra));
        p=voigt(x,pos,wid,extra);
    case 35
        p=GL(x,pos,wid,extra);
    case 36
        p=expgaussian(x,pos,wid,-extra);
    case 37
        p=pearson(x,pos,wid,extra);
    case 38
        p=explorentzian(x,pos,wid,-extra);
    case 40
        p=sine(x,pos,wid);
    case 41
        p=rectangle(x,pos,wid);
    case 42
        p=ngaussian(x,pos,wid,extra);
    case 43
        p=Gompertz(x,pos,wid,extra);
    case 44
        p=OneMinusExp(x,pos,wid);
    case 45
        p=FourPL(x,pos,wid,extra);
    case 46
        p=quadslope(x,pos,wid);
    case 47
        p=blackbody(x,pos);
    case 48
        p=exppulse(x,pos,wid);
    otherwise
end % switch
% ----------------------------------------------------------------------
function [coef,RSquared,BootResults]=plotit(xi,yi,polyorder)
% plotit accepts data in the form of a single vector, a pair of vectors,
% "x" and "y", or a 2xn or nx2 matrix with x in first row or column. If the
% input argument "polyorder" is supplied, it fits the data to a polynomial
% of order "polyorder", plots the data in red dots and the fit as a blue
% line, and displays the fit coefficients and R-squared in the upper left
% corner of the graph. Polyorder=1 for straight line, =2 for quadratic
% (parabola) etc. If the 3rd output argument (BootResults) is supplied,
% computes coefficient error estimates by the bootstrap method and returns
% the results in the matrix "BootResults" (of size 5 x polyorder+1). You
% can change the number of bootstrap samples in line 48. 
%
% Copyright (c) 2012, Thomas C. O'Haver

% process input arguments
NumTrials=100;
if nargin==1, % Single argument, might be y vector or x,y matrix
    polyorder=0;
    datasize=size(xi);
    if isvector(xi), % xi is vector, not matrix
        y=xi;
        x=1:length(y); % Use this only to create an x vector if needed
        xlabelstring='index number';
    else % x is a matrix, split it into x and y vectors
        if datasize(1)<datasize(2),xi=xi';end
        x=xi(:,1);
        y=xi(:,2);
        polyorder=0;
        xlabelstring='x';
    end
end

if nargin==2, % Two arguments, might be x,y vectors or matrix and polyorder
    datasize=size(xi);
    if isscalar(yi),   %  matrix and polyorder
        if datasize(1)<datasize(2),xi=xi';end
        if isvector(xi),
            x=1:length(xi);
            xlabelstring='index number';
            y=xi';
            polyorder=yi;
        else
            x=xi(:,1);
            y=xi(:,2);
            polyorder=yi;
            xlabelstring='x';
        end
    else %  x,y vectors
        x=xi;
        y=yi;
        polyorder=0;
        xlabelstring='x';
    end
end

if nargin==3, % Tree arguments, must be x,y vectors and polyorder
    x=xi;
    y=yi;
    xlabelstring='x';
end

x=reshape(x,[1,length(x)]);
y=reshape(y,[1,length(y)]);

if length(x)==length(y),
    if polyorder>0, % Plot and fit the data
        % Compute the fit
        coef=polyfit(x,y,polyorder);
        % Plot the data and the polynomial fit
        xx=linspace(min(x),max(x));
        yhat=polyval(coef,xx);
        plot(x,y,'.b',xx,yhat,'-r')
        axis([min(x) max(x) min(y) max(y)]);
        xlabel(xlabelstring);ylabel('y')
        title(['Number of data points= ' num2str(length(x)) ] )
        % Compute the correlation coefficient and R-Squared
        if IsOctave,
            cc=corr(polyval(coef,x'),y');
            RSquared=cc.^2;
        else
            cc=corrcoef(polyval(coef,x),y);
            RSquared=cc(2).^2;
        end %   if IsOctave,
        % Label the graph with the fit information
        text(min(x),max(y)-.1.*(max(y)-min(y)),['   Polynomial Order of fit = ' num2str(polyorder)] );
        text(min(x),max(y)-.2.*(max(y)-min(y)),['   Fit coefficients = ' num2str(coef)] );
        text(min(x),max(y)-.3.*(max(y)-min(y)),['   R-Squared = ' num2str(RSquared)] );
    else % Just plot the data, no fit.
        plot(x,y)
        axis([min(x) max(x) min(y) max(y)]);
        xlabel(xlabelstring);ylabel('y')
        title(['Number of data points= ' num2str(length(x)) ] )
        coef=[];
        RSquared=[];
    end
else
    disp('Error: x and y must be the same size.')
    sizex=size(x)
    sizey=size(y)
end
    
if nargout==3, % Compute the bootsrtap
    BootstrapResultsMatrix=zeros(polyorder+1,NumTrials);
    clear bx by
    for trial=1:NumTrials,
        n=1;
        bx=x;
        by=y;
        while n<length(x)-1,
            if rand>.5,
                bx(n)=x(n+1);
                by(n)=y(n+1);
            end % if rand>.5,
            n=n+1;
        end % while n<length(xx)-1,
        coef=polyfit(bx,by,polyorder);
        BootstrapResultsMatrix(:,trial)=coef;
    end % for trial=1:NumTrials,
    disp(' ')
    disp('Bootstrap Results');
    BootstrapMean=mean(real(BootstrapResultsMatrix(:,:)'));
    BootstrapSTD=std(BootstrapResultsMatrix(:,:)');
    BootstrapIQR=iqr(BootstrapResultsMatrix(:,:)');
    PercentRSD=100.*BootstrapSTD./abs(BootstrapMean);
    PercentIQR=100.*BootstrapIQR./abs(BootstrapMean);
    BootstrapMean=BootstrapMean(1:polyorder+1);
    BootstrapSTD=BootstrapSTD(1:polyorder+1);
    BootstrapIQR=BootstrapIQR(1:polyorder+1);
    PercentRSD=PercentRSD(1:polyorder+1);
    PercentIQR=PercentIQR(1:polyorder+1);
        disp(['Mean:       ', num2str(BootstrapMean)])
        disp(['STD:        ', num2str(BootstrapSTD)])
        disp(['RSD (IQR):  ', num2str(BootstrapIQR)])
        disp(['% RSD:      ', num2str(PercentRSD)])
        disp(['% RSD(IQR): ', num2str(PercentIQR) '%' ])
    BootResults(1,:)=BootstrapMean;
    BootResults(2,:)= BootstrapSTD ;
    BootResults(3,:)=PercentRSD;
    BootResults(4,:)=BootstrapIQR;
    BootResults(5,:)=PercentIQR;
end % if NumArgOut==3,

function isOctave = IsOctave()
% Returns true if this code is being executed by Octave.
% Returns false if this code is being executed by MATLAB, or any other MATLAB
% variant.
% 
%    usage: isOctave = IsOctave()
    persistent octaveVersionIsBuiltIn;
    if (isempty(octaveVersionIsBuiltIn))
        octaveVersionIsBuiltIn = (exist('OCTAVE_VERSION', 'builtin') == 5);
        % exist returns 5 to indicate a built-in function.
    end
    isOctave = octaveVersionIsBuiltIn;
    % If OCTAVE_VERSION is a built-in function, then we must be in Octave.
    % Since the result cannot change between function calls, it is cached in a
    % persistent variable.  isOctave cannot be a persistent variable, because it
    % is the return value of the function, so instead the persistent result must
    % be cached in a separate variable
    
