function ipf(arg1,arg2,arg3,arg4)
% Keyboard-operated Interactive Peak Fitter for data in separate x,y
% vectors, or in a single y vector, or in a 2xn or nx2 data matrix with
% x values in row/column 1 and y values in row/column 2 (e.g. [x y]).
% Press K key for list of keyboard commands.
% Version 10.3, May, 2014, replaces sigmoid (logistic function) with
% separate "up sigmoid" (shape 10) and "down sigmoid" (shape 23).
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
% See http://terpconnect.umd.edu/~toh/spectrum/CurveFittingC.html
% See http://www.wam.umd.edu/~toh/spectrum/InteractivePeakFitter.htm
% T. C. O'Haver (toh@umd.edu).
%
% Example 1: 
% x=[0:.005:1];y=humps(x).^3;ipf(x,y)
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
% Select peak shape from menu - (minus or hyphen), then type number and Enter
% Select peak shape...........g Gaussian
%                             G (Shift-G) Fixed-width Gaussian
%                             Shift-P  Fixed-position Gaussians
%                             h Equal-width Gaussian
%                             H (Shift-H) bifurcated Gaussian (a,z keys adjust shape)
%                             l Lorentzian
%                             ; Equal-width Lorentzian
%                             L (Shift-L) Fixed-width Lorentzian
%                             Shift-[  Fixed-position Lorentzians
%                             : bifurcated Lorentzians (a,z keys adjust shape)
%                             o Logistic distribution (See "S" for logistic function)
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
% Fit.........................f
% Select autozero mode........t  selects none, linear, quadratic, or flat baseline mode
% Monopolar/bipolar mode......+  Flips between + peaks only and +/- peak mode
% Toggle log y mode OFF/ON....m  Plot linear or log Y axis on lower graph
% 2-point Baseline............b, then click left and right baseline
% Set manual baseline.........Backspace, then click baseline at multiple points
% Restore original baseline...\  to cancel previous background subtraction
% Cursor start position.......c, then click on each peak position
% Type in start vector........C (Shift-C) Type or Paste start vector [p1 w1 p2 w2 ...]
% Enter value of 'Extra'......Shift-x, type value, press Enter.
% Adjust 'Extra' up/down......a,z`: 5% change; upper case A,Z: 0.5% change.
% Print parameters & results..q
% Print fit results only......r
% Plot signal in figure 2.....y
% Print out model data table..d
% Refine fit..................x  Takes best of 10 trial fits (change in line 220)
% Print peakfit function......w  with all input arguments
% Compute bootstrap stats.....v  Estimates standard deViations of parameters.
% Fit single bootstrap........n  Fits siNgle bootstrap sub-sample
% Save Figure as png file.....Shift-S  Saves as Figure1.png, Figure2.png....
% Display current settings....Shift-? displays table of current values
% Prints list of commands.....k

% Copyright (c) 2013, Thomas C. O'Haver
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
extra=length(X)/100.*zeros(size(Shape));
if Shape==20,extra=1;end
FigNum=1;
delta=0; % Initial Random change in initial start guesses for each re-fit
% xo=length(Y)/2; % Initial Pan setting
% dx=length(Y)/4; % Initial Zoom setting
PEAKHEIGHTS=zeros(NumPeaks,1);
xxx=zeros(1,200);
AA=zeros(NumPeaks,200);
AUTOZERO=0; % Default to no baseline correction Press T to select mode.
logplot=0; % Default to linear mode. Press m to chage to log.
SavedSignal=Y;
FIXEDPARAMETERS=0; % Used only for fixed position or width shapes.
BIPOLAR=0; % Default to monopolar (positive peaks only) mode.

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
global  FIXEDPARAMETERS logplot et FigNum BIPOLAR shapesvector
% (et=elapsed trie)
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
            if logplot==0,logplot=1;else logplot=0;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 27 % When 'ESC' key is pressed, resets pan and zoom
            xo=length(Y)/2; % Initial Pan setting
            dx=length(Y)/4; % Initial Zoom setting
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 102
            % When 'f' key is pressed, tweaks start values, computes and
            % plots fit.
            startnow=start;
            delta=(max(xx)-min(xx))/100;
            for k=1:2*NumPeaks,
                startnow(k)=start(k)+(rand-.5)*delta;
            end
            [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,startnow,extra);
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
                peakfit([xx;yy],0,0,NumPeaks,shapesvector,extra,1,start,AUTOZERO,FIXEDPARAMETERS,1);
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
        case {103,108,111,112,101,104,59,106,117,115,71,76,96,72,58,80,123,69,85,86,68}
            % Selects peak shape when the following keys are pressed.
            switch key
                case 103 % When 'g' key is pressed, peak shape is set to Gaussian.
                    n=1;
                case 108 % When 'l' key is pressed, peak shape is set to Lorentzian.
                    n=2;
                case 111 % When 'o' key is pressed, peak shape is set to Logistic distribution.
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
                case 58 % When ':' key is pressed, peak shape is set to bifurcated Lorentzian
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
            disp('Gaussian: y=exp(-((x-pos)./(0.6005615.*width)) .^2)')
            disp('   Gaussians with independent positions and widths : 1 (default)')
            disp('   Exponentional-broadened Gaussian : 5 ')
            disp('   Exponentional-broadened equal-width Gaussian : 8')
            disp('   Gaussians with the same widths : 6')
            disp('   Gaussians with preset fixed widths : 11')
            disp('   Fixed-position Gaussians : 16 ')
            disp('   Asymmetrical Gaussians with unequal half-widths on both sides : 14')      
            disp('Lorentzian: y=ones(size(x))./(1+((x-pos)./(0.5.*width)).^2)')
            disp('   Lorentzians with independent positions and widths : 2')
            disp('   Exponentional-broadened Lorentzian : 18 ')            
            disp('   Equal-width Lorentzians : 7  ')
            disp('   Fixed-width Lorentzian : 12')
            disp('   Fixed-position Lorentzian : 17')
            disp('   Asymmetrical Lorentzians with unequal half-widths on both sides : 15')
            disp('Blended sum of Gaussian and Lorentzian functions : 13')
            disp('Voigt profile: 20')            
            disp('Logistic distribution: n=exp(-((x-pos)/(.477.*wid)).^2); y=(2.*n)./(1+n) : 3  ')
            disp('Pearson: y=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m : 4')
            disp('Exponential pulse: y=(x-tau2)./tau1.*exp(1-(x-tau2)./tau1) : 9')
            disp('Alpha function: y=(x-spoint)./pos.*exp(1-(x-spoint)./pos); : 19')           
            disp('Up Sigmoid (logistic function): y=.5+.5*erf((x-tau1)/sqrt(2*tau2)) : 10')
            disp('Down Sigmoid y=.5-.5*erf((x-tau1)/sqrt(2*tau2) ): 23')
            disp('For multi-shape, type numbers in brackets: e.g. [1 1 2 4] ')
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
            if Shape>22, Shape=22;end
            if Shape<1, Shape=1;end
            if isempty(Shape),Shape=1;end
            switch Shape
                case 1
                    ShapeString='Gaussian';
                case 2
                    ShapeString='Lorentzian';
                case 3
                    ShapeString='Logistic distribution';
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
                    inputwidth=input('Peak width: ');
                    if isempty(inputwidth),
                    else
                        FIXEDPARAMETERS=inputwidth;
                    end
                case 12
                    ShapeString='Fixed-width Lorentzian';
                    inputwidth=input('Peak width: ');
                    if isempty(inputwidth),
                    else
                        FIXEDPARAMETERS=inputwidth;
                    end
                case 13
                    ShapeString='Gauss/Lorentz blend';
                case 14
                    ShapeString='bifurcated Gaussian';
                case 15
                    ShapeString='bifurcated Lorentzian';
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
                otherwise
            end % switch
            figure(1);subplot(2,1,2)
            if BIPOLAR,
                xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = + - ' ShapeString  ] )
            else
                xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = + ' ShapeString  ] )
            end
        case 88 % Shift-X
            disp(['Current calue of extra parameter: ' num2str(extra) ] )
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
            disp(['Fitting Error = ' num2str(MeanFitError) '%'  ])
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
        case 116
            % When 't' key is pressed, steps through AUTOZERO modes
            AUTOZERO=AUTOZERO+1;
            if AUTOZERO==4,AUTOZERO=0;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 107
            % When 'k' key is pressed, prints out table of keyboard commands
            disp('KEYBOARD CONTROLS (Version 10.2):')
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
            disp('                             G (Shift-G) fixed-width Gaussians')
            disp('                             Shift P fixed-position Gaussians')
            disp('                             H (Shift-H) bifurcated Gaussians (a,z keys adjust shape)')
            disp('                             e Exponential-broadened Gaussian')
            disp('                             j exponential-broadened equal-width Gaussians')
            disp('                                 (a,z keys adjust broadening)')
            disp('                             l Lorentzian')
            disp('                             ; equal-width Lorentzians')
            disp('                             Shift [ fixed-position Lorentzians')
            disp('                             : bifurcated Lorentzians (a,z keys adjust shape)')
            disp('                             E (Shift-E) Exponential-broadened Lorentzians')
            disp('                             L (Shift-L) Fixed-width Lorentzians')
            disp('                             o Logistic distribution')
            disp('                             p Pearson (a,z keys adjust shape)')
            disp('                             u exponential pUlse  y=exp(-tau1.*x).*(1-exp(-tau2.*x))')
            disp('                             U (Shift-U) Alpha: y=(x-tau2)./tau1.*exp(1-(x-tau2)./tau1)')
            disp('                             Up Sigmoid (logistic function): y=.5+.5*erf((x-tau1)/sqrt(2*tau2))')
            disp('                             Shift-D Down Sigmoid y=.5-.5*erf((x-tau1)/sqrt(2*tau2))')
            disp('                             ~ Gauss/Lorentz blend (a/z keys adjust shape)')
            disp('                             V (Shift-V) Voigt profile (a/z adjusts shape)')
            disp(' Fit.........................f Performs one fit of the specified model')
            disp(' Select autozero mode........t  selects none, linear, quadratic, or flat baseline mode')
            disp(' Monopolar/bipolar mode......+  Flips between + peaks only and +/- peak mode')
            disp(' Toggle log y mode OFF/ON....m  Plot linear or log Y axis on lower graph')
            disp(' 2-point Baseline............b, then click left and right baseline')
            disp(' Set manual baseline.........Backspace, then click baseline at multiple points')
            disp(' Restore original baseline...\  to cancel previous background subtraction')
            disp(' Cursor start positions......c, then click on each peak position')
            disp(' Type in start vector........C (Shift-C) Type or Paste start vector [p1 w1 p2 w2 ...]')
            disp(' Enter value of ''Extra''......Shift-x, type value or [vector], press Enter.')
            disp(' Adjust ''Extra'' up/down......a,z: 5% change; upper case A,Z: 0.5% change.')
            disp(' Print parameters & results..q')
            disp(' Print fit results only......r')
            disp(' Compute bootstrap stats.....v  Estimates standard deViations of parameters.')
            disp(' Fit single bootstrap........n  Extracts and Fits siNgle bootstrap sub-sample.')
            disp(' Plot signal in figure 2.....y')
            disp(' Print model data table......d')
            disp(' Refine fit..................x Takes best of 10 trial fits (change in line 177)')
            disp(' Print peakfit function......w  Print peakfit function with all input arguments')
            disp(' Save Figure as png file.....Shift-S  Saves as Figure1.png, Figure2.png, etc.')
            disp(' Display current settings....Shift-? displays table of current values')
        case 120 % When 'x' key is pressed, calls peakfit to take best of 'NumTrials' trial fits
            center=(max(xx)+min(xx))/2;
            window=max(xx)-min(xx);
            disp(['Calculating best of ' num2str(NumTrials) ' trial fits....' ])
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
            [FitResults,MeanFitError]=peakfit([X',Y'],center,window,NumPeaks,shapesvector,extra,NumTrials,FirstGuess,AUTOZERO,fixedparameters);
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
                disp(['ipf(datamatrix,' num2str(center) ',' num2str(window) ')']);
                if Shape==11||12||16||17,
                    disp(['[FitResults,MeanFitError]=peakfit(datamatrix,' num2str(center) ',' num2str(window) ',' num2str(NumPeaks) ',' num2str(Shape) ',' num2str(extra) ',' num2str(NumTrials) ', [' num2str(FirstGuess) '], ' num2str(AUTOZERO) ', [' num2str(FIXEDPARAMETERS) '] ) ' ]);
                else
                    disp(['[FitResults,MeanFitError]=peakfit(datamatrix,' num2str(center) ',' num2str(window) ',' num2str(NumPeaks) ',' num2str(Shape) ',' num2str(extra) ',' num2str(NumTrials) ', [' num2str(FirstGuess) '], ' num2str(AUTOZERO) ' ) ']);
                end
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
            if BIPOLAR,
                disp('Bipolar mode')
            else
                disp('Positive peaks only')
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
            switch Shape
                case 4
                    disp(['Shape Constant = ' num2str(extra)])  
                case {5,8,18}
                    disp(['Time Constant = ' num2str(extra)])
                case 13
                    disp(['% Gaussian = ' num2str(extra)])
                case {14,15}
                    disp(['Asymmetry = ' num2str(extra)])
                case 20
                    disp(['Alpha = ' num2str(extra)])    
            end
            if Shape==11||Shape==12, disp(['Fixed Peak Width = ' num2str(FIXEDPARAMETERS)]), end
            disp(['Fitted x range = ' num2str(min(xx)) ' - ' num2str(max(xx)) ' (dx=' num2str(max(xx)-min(xx)) ')  (Center=' num2str((max(xx)+min(xx))/2) ')  ' ])
            disp([num2str(length(xx)) ' data points fit' ])
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
        case 121    % When 'Y' key is pressed (Added on version 5)
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
            disp(['Estimated time: ',num2str(EstTime) ' seconds']);
            if NumTrialsBoot,
                disp('Computing bootstrap sampling statistics....May take several minutes.')
                tic;
                BootstrapResultsMatrix=zeros(5,100,NumPeaks);
                BootstrapErrorMatrix=zeros(1,100,NumPeaks);
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
                    [FitResults,BootFitError]=peakfit([bx,by],center,window,NumPeaks,Shape,extra,NumTrialsBoot,start,AUTOZERO,FIXEDPARAMETERS,1);
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
                    disp(['Bootstrap Mean: ', num2str(BootstrapMean(2:5))])
                    disp(['Bootstrap STD:  ', num2str(BootstrapSTD(2:5))])
                    disp(['Bootstrap IQR:  ', num2str(BootstrapIQR(2:5))])
                    disp(['Percent RSD:    ', num2str(PercentRSD(2:5))])
                    disp(['Percent IQR:    ', num2str(PercentIQR(2:5))])
                    % mean(PercentIQR(2:5)./PercentRSD(2:5))
                end
                toc;
                disp(['Min/Max Fitting Error of the 100 bootstrap samples: ', num2str(MaxError) '/' num2str(MinError) ])
                disp('-------------------------------------------------------------------')
            end
            figure(1)
            title('ipf 10  Typical Bootstrap sample fit')
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
            [FitResults,MeanFitError]=peakfit([bx,by],center,window,NumPeaks,shapesvector,extra,NumTrialsBoot,start,AUTOZERO,FIXEDPARAMETERS,1);
%             disp(['Fitting Error = ' num2str(MeanFitError) '%'])
%             if Shape==9||Shape==10||Shape==19,
%                 disp('         Peak#     Tau1         Height       Tau2          Area');
%             else
%                 disp('         Peak#   Position       Height       Width         Area');
%             end
%             disp(FitResults)
            figure(1)
            title('ipf 10  Single Bootstrap sample fit')
        case 83
            saveas(gcf,['Figure' num2str(FigNum) '.png'])
            FigNum=FigNum+1;
        case 63 % Shift-? Displays current settings
            disp(' ')
            disp( 'Current Settings:')
            disp([ 'Shape = ' num2str(Shape) ])
            disp([ 'BIPOLAR = ' num2str(BIPOLAR) ])
            disp([ 'shapesvector = ' num2str(shapesvector) ])
            disp([ 'NumPeaks = ' num2str(NumPeaks) ])
            disp([ 'Baseline mode = ' num2str(AUTOZERO) ])
            disp(['Fitted x range = ' num2str(min(xx)) ' - ' num2str(max(xx)) ' (dx=' num2str(max(xx)-min(xx)) ')  (Center=' num2str((max(xx)+min(xx))/2) ')  ' ])
            disp([ 'start vector = [' num2str(start) ']' ])
            disp([ 'extra = ' num2str(extra) ])
            disp([ 'NumTrials = ' num2str(NumTrials) ])
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
global AUTOZERO AA xxx PEAKHEIGHTS logplot BIPOLAR
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
plot(xx,yy,'b.'); % Plot the original signal in blue in upper window
xlabel('Line up the estimated peak positions roughly with the vertical lines')
lyy=min(yy);
uyy=max(yy)+(max(yy)-min(yy))/10;
switch AUTOZERO,
    case 0
        title('ipf 10 No baseline correction. Pan and Zoom to isolate peaks to be fit in upper window.')
        if lyy<uyy;axis([X(Startx) X(Endx) lyy uyy ]);end
    case 1
        title('ipf 10 Linear baseline mode. Pan and Zoom to isolate peaks to be fit in upper window.')
        if lyy<uyy;axis([X(Startx) X(Endx) lyy uyy ]);end
    case 2
        title('ipf 10 Quadratic baseline mode. Pan and Zoom to isolate peaks to be fit in upper window.')
        if lyy<uyy;axis([X(Startx) X(Endx) lyy uyy ]);end
    case 3
        title('ipf 10 Flat baseline mode. Pan and Zoom to isolate peaks to be fit in upper window.')
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
    plot([markx markx],[lyy uyy],'m--')
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
global PEAKHEIGHTS AUTOZERO AA xxx residual FIXEDPARAMETERS BIPOLAR 
global shapesvector
PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
x0=min(xx);
% StartVector=start
% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.00001,'TolFun',.00001,'Display','off','MaxFunEvals',1000);
tic
switch Shape
    case 1
        FitParameters=fminsearch(@fitgaussian,start,options,xx,yy);
        ShapeString='Gaussian';
    case 2
        FitParameters=fminsearch(@fitlorentzian,start,options,xx,yy);
        ShapeString='Lorentzian';
    case 3
        FitParameters=fminsearch(@fitlogistic,start,options,xx,yy);
        ShapeString='Logistic';
    case 4
        FitParameters=fminsearch(@fitpearson,start,options,xx,yy,extra(1));
        ShapeString='Pearson';
    case 5
        FitParameters=fminsearch(@fitexpgaussian,start,options,xx,yy,-extra(1));
        ShapeString='ExpGaussian';
    case 6
        cwnewstart(1)=start(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=start(2.*pc-1);
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        FitParameters=fminsearch(@fitewgaussian,cwnewstart,options,xx,yy);
        ShapeString='Equal width Gauss.';
    case 7
        cwnewstart(1)=start(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=start(2.*pc-1);
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        FitParameters=fminsearch(@fitewlorentzian,cwnewstart,options,xx,yy);
        ShapeString='Equal width Lorentzian';
    case 8
        cwnewstart(1)=start(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=start(2.*pc-1);
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        FitParameters=fminsearch(@fitexpewgaussian,cwnewstart,options,xx,yy,-extra(1));
        ShapeString='Exp. equal width Gaussians';
    case 9
        FitParameters=fminsearch(@fitexppulse,start,options,xx,yy);
        ShapeString='Exponential Pulse';
    case 10
        FitParameters=fminsearch(@fitupsigmoid,start,options,xx,yy);
        ShapeString='Up Sigmoid (logistic function)';
    case 23
        FitParameters=fminsearch(@fitdownsigmoid,start,options,xx,yy);
        ShapeString='Down Sigmoid (logistic function)';
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
        FitParameters=fminsearch(@FitGL,start,options,xx,yy,extra(1));
        ShapeString='Gausss/Lorentz blend';
    case 14
        FitParameters=fminsearch(@fitBiGaussian,start,options,xx,yy,extra(1));
        ShapeString='bifurcated Gaussian';
    case 15
        FitParameters=fminsearch(@fitBiLorentzian,start,options,xx,yy,extra(1));
        ShapeString='Alpha function';
    case 16
        fixedstart=[];
        for pc=1:NumPeaks,
            fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
            fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc);
        end
        FitParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart,options);
        ShapeString='Fixed-position Gaussians';
    case 17
        fixedstart=[];
        for pc=1:NumPeaks,
            fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
            fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc);
        end
        FitParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart,options);
        ShapeString='Fixed-position Lorentzians';
    case 18
        FitParameters=fminsearch(@fitexplorentzian,start,options,xx,yy,-extra(1));
        ShapeString='ExpLorentzian';
    case 19
        FitParameters=fminsearch(@fitalphafunction,start,options,xx,yy);
        ShapeString='Alpha function';
    case 20
        FitParameters=fminsearch(@fitvoigt,start,options,xx,yy,extra(1));
        ShapeString='Voigt profile';
    case 21
        FitParameters=fminsearch(@(lambda)(fittriangular(lambda,xx,yy)),start,options);
        ShapeString='Triangular';
    case 22
        FitParameters=fminsearch(@(lambda)(fitmultiple(lambda,xx,yy,shapesvector,extra)),start,options);
        ShapeString=num2str(shapesvector);
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
            A(m,:)=logistic(xx,FitParameters(2*m-1),FitParameters(2*m));
            AA(m,:)=logistic(xxx,FitParameters(2*m-1),FitParameters(2*m));
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
            A(m,:)=BiLorentzian(xx,FitParameters(2*m-1),FitParameters(2*m),extra(1));
            AA(m,:)=BiLorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra(1));
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
figure(1);subplot(2,1,1);plot(xx,yy,'b.'); % Plot the original signal in blue dots
hold on
for m=1:NumPeaks,
    plot(xxx,Heights(m)*AA(m,:),'g')  % Plot the individual component peaks in green lines
    area(m)=trapz(xx,Heights(m)*A(m,:)); % Compute the area of each component peak using trapezoidal method
end
% Mark starting peak positions with vertical dashed lines
for marker=1:NumPeaks,
    markx=start((2*marker)-1);
    subplot(2,1,1);plot([markx markx],[0 max(yy)],'m--')
end % for
plot(xxx,mmodel,'r');  % Plot the total model (sum of component peaks) in red lines
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
        title('ipf 10 No baseline correction. Pan and Zoom to isolate peaks to be fit in upper window.')
    case 1
        title('ipf 10 Linear baseline mode. Pan and Zoom to isolate peaks to be fit in upper window.')
    case 2
        title('ipf 10 Quadratic baseline mode. Pan and Zoom to isolate peaks to be fit in upper window.')
    case 3
        title('ipf 10 Flat baseline mode. Pan and Zoom to isolate peaks to be fit in upper window.')
end
xlabel('Vertical dotted lines indicate first guess peak positions. C to customize.');
% Bottom half of the figure shows the residuals and displays RMS error
% between original signal and model
residual=yy-model;
figure(1);subplot(2,1,2);plot(xx,residual,'b.')
MeanFitError=100*norm(residual)./(sqrt(n)*max(yy));
title('F: fit again, X: best of 10 trials, T: baseline, Shift-X: enter Extra, a/z adjust Extra, Q: report, R: peak table')
ylabel('Residuals')
switch Shape
    case 4
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '   Shape Constant = ' num2str(extra) '   Error = ' num2str(round(1000*MeanFitError)/1000) '%' ] )
    case {5,8,18}
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '    Time Constant = ' num2str(extra(1)) '   Error = ' num2str(round(1000*MeanFitError)/1000) '%' ] )
    case 13
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '       % Gaussian = ' num2str(extra(1)) '   Error = ' num2str(round(1000*MeanFitError)/1000) '%' ] )
    case {14,15}
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '        Asymmetry = ' num2str(extra(1)) '   Error = ' num2str(round(1000*MeanFitError)/1000) '%' ] )
    case 20
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '            Alpha = ' num2str(extra(1)) '   Error = ' num2str(round(1000*MeanFitError)/1000) '%' ] )
    otherwise
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '        Error = ' num2str(round(1000*MeanFitError)/1000) '% ' ] )
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

% Display Fit Results on upper graph
subplot(2,1,1);
startx=min(xx)+(max(xx)-min(xx))./20;
dxx=(max(xx)-min(xx))./10;
dyy=(max(yy)-min(yy))./10;
starty=max(yy)-dyy;
FigureSize=get(gcf,'Position');
if Shape==9||Shape==10, % pulse and sigmoid shapes
    text(startx,starty+dyy/2,['Peak #          tau1           Height           tau2             Area'] );
else
    if Shape==11||Shape==12, % Fixed-width shapes
        text(startx,starty+dyy/2,['Peak #          Position        Height         width (fixed)    Area'] );
        
    else
        text(startx,starty+dyy/2,['Peak #          Position        Height         Width             Area'] );
    end
end
% Alternative method of printing FitResults using sprintf
for peaknumber=1:NumPeaks,
    if Shape==0,LastColumn=4; else LastColumn=5;end
    for column=1:LastColumn,
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
% Alternative method of printing FitResults using num2str
% text(startx, starty-(NumPeaks./2).*dyy,num2str(FitResults));
% ----------------------------------------------------------------------
function [FitResults,LowestError,baseline,BestStart,xi,yi,BootResults]=peakfit(signal,center,window,NumPeaks,Shape,extra,NumTrials,start,autozero,fixedparameters,plots)
global AA xxx PEAKHEIGHTS FIXEDWIDTH FIXEDPOSITIONS AUTOZERO delta c BIPOLAR
% peakfit.m modified versison of version 5.1, February 2014
format short g
format compact
warning off all
plots=1;
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

% Default values for placeholder zeros
peakshape=Shape;
if NumTrials==0;NumTrials=1;end
if isscalar(peakshape),
else
    % disp('peakshape is vector');
    shapesvector=peakshape;
    NumPeaks=length(shapesvector);
    peakshape=22;
end
if peakshape==0;peakshape=1;end
if NumPeaks==0;NumPeaks=1;end
if start==0;start=calcstart(xx,NumPeaks,xoffset);end
if FIXEDWIDTH==0, FIXEDWIDTH=length(xx)/10;end
if peakshape==16;FIXEDPOSITIONS=fixedparameters;end
if AUTOZERO>3,AUTOZERO=3,end
if AUTOZERO<0,AUTOZERO=0,end
delta=1;

% % Remove linear baseline from data segment if AUTOZERO==1
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
        ShapeString='Logistic distribution';
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
        ShapeString='BiLorentzian';   
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
    % StartVector=newstart
    % disp(['Trial number ' num2str(k) ] ) % optionally prints the current trial number as progress indicator
  switch peakshape(1)
    case 1
        TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart,options);
    case 2
        TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart,options);
    case 3
        TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart,options);
    case 4
        TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart,options);
    case 5
        zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
        zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
        TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart,options);
    case 6
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart,options);
    case 7
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitlorentziancw(lambda,xx,yy)),cwnewstart,options);
    case 8
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart,options);
      case 9
          TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart,options);
      case 10
          TrialParameters=fminsearch(@(lambda)(fitupsigmoid(lambda,xx,yy)),newstart,options);
      case 23
          TrialParameters=fminsearch(@(lambda)(fitdownsigmoid(lambda,xx,yy)),newstart,options);
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
      case 14
          TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart,options);
      case 15
          TrialParameters=fminsearch(@(lambda)(fitBiLorentzian(lambda,xx,yy,extra)),newstart,options);
      case 16
           fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
              fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart,options);
      case 17
           fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
              fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart,options);
      case 18
        zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
        zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ];
%         plot(zxx,zyy);pause
        TrialParameters=fminsearch(@(lambda)(fitexplorentzian(lambda,zxx,zyy,-extra)),newstart,options);     
      case 19
          TrialParameters=fminsearch(@(lambda)(fitalphafunction(lambda,xx,yy)),newstart,options);
      case 20
          TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart,options);
      case 21
          TrialParameters=fminsearch(@(lambda)(fittriangular(lambda,xx,yy)),newstart,options);
      case 22
          TrialParameters=fminsearch(@(lambda)(fitmultiple(lambda,xx,yy,shapesvector,extra)),newstart,options);
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
            A(m,:)=gaussian(xx,TrialParameters(m),FIXEDWIDTH);
        case 12
            A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDWIDTH);
        case 13
            A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 14
            A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 15
            A(m,:)=BiLorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
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
        otherwise
    end % switch
    
    for parameter=1:2:2*NumPeaks,
        newstart(parameter)=newstart(parameter)*(1+delta*(rand-.5)/500);
        newstart(parameter+1)=newstart(parameter+1)*(1+delta*(rand-.5)/100);
    end
end % for NumPeaks

% Multiplies each row by the corresponding amplitude and adds them up
if AUTOZERO==3,
    baseline=PEAKHEIGHTS(1);
    Heights=PEAKHEIGHTS(2:1+NumPeaks);
    model=Heights'*A+baseline;
else
    model=PEAKHEIGHTS'*A;
    Heights=PEAKHEIGHTS;
    baseline=0;
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
end % for k (NumTrials)
%
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
        AA(m,:)=gaussian(xxx,FitParameters(m),FIXEDWIDTH);
    case 12
        AA(m,:)=lorentzian(xxx,FitParameters(m),FIXEDWIDTH);
    case 13
        AA(m,:)=GL(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 14
        AA(m,:)=BiGaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 15
        AA(m,:)=BiLorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
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
       otherwise
  end % switch
end % for NumPeaks

% Multiplies each row by the corresponding amplitude and adds them up
heightsize=size(height');
AAsize=size(AA);
if heightsize(2)==AAsize(1),
   mmodel=height'*AA+baseline;
else
    mmodel=height*AA+baseline;
end
% Top half of the figure shows original signal and the fitted model.
if plots,
    subplot(2,1,1);plot(xx+xoffset,yy,'b.'); % Plot the original signal in blue dots
    hold on
end
for m=1:NumPeaks,
    if plots, plot(xxx+xoffset,height(m)*AA(m,:)+baseline,'g'),end  % Plot the individual component peaks in green lines
    area(m)=trapz(xxx+xoffset,height(m)*AA(m,:)); % Compute the area of each component peak using trapezoidal method
    yi(m,:)=height(m)*AA(m,:); % (NEW) Place y values of individual model peaks into matrix yi
end
xi=xxx+xoffset; % (NEW) Place the x-values of the individual model peaks into xi

if plots,
    % Mark starting peak positions with vertical dashed lines
    if peakshape(1)==16||peakshape(1)==17
    else
        for marker=1:NumPeaks,
            markx=BestStart((2*marker)-1);
            subplot(2,1,1);plot([markx+xoffset markx+xoffset],[0 max(yy)],'m--')
        end % for
    end % if peakshape
    plot(xxx+xoffset,mmodel,'r');  % Plot the total model (sum of component peaks) in red lines
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
            title(['peakfit 5   No baseline correction'])
        case 1
            title(['peakfit 5   Linear baseline subtraction'])
        case 2
            title(['peakfit 5   Quadratic subtraction baseline'])
        case 3
            title(['peakfit 5   Flat baseline correction'])
    end
 
    switch peakshape(1)
    case {4,20}
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '      Shape Constant = ' num2str(extra) '   Error = ' num2str(round(1000*LowestError)/1000) '%' ] )
    case {5,8,18}
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '      Time Constant = ' num2str(extra) '   Error = ' num2str(round(1000*LowestError)/1000) '%' ] )
    case 13
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '      % Gaussian = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '% ' ] )
    case {14,15,22}
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '      extra = ' num2str(extra) '     Error = ' num2str(round(1000*LowestError)/1000)  '% ' ] )
    otherwise
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(1000*LowestError)/1000) '% ' ] )
    end

    % Bottom half of the figure shows the residuals and displays RMS error
    % between original signal and model
    residual=yy-bestmodel;
    subplot(2,1,2);plot(xx+xoffset,residual,'b.')
    axis([min(xx)+xoffset max(xx)+xoffset min(residual) max(residual)]);
    xlabel('Residual Plot')
end % if plots

% Put results into a matrix, one row for each peak, showing peak index number,
% position, amplitude, and width.
for m=1:NumPeaks,
    if m==1,
        if peakshape(1)==6||peakshape(1)==7||peakshape(1)==8, % equal-width peak models only
            FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape(1)==11||peakshape(1)==12, % Fixed-width shapes only
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                if peakshape(1)==16||peakshape(1)==17, % Fixed-position shapes only
                    FitResults=[round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)];
                else
                    FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)];
                end
            end
        end % if peakshape
    else
        if peakshape(1)==6||peakshape(1)==7||peakshape(1)==8, % equal-width peak models only
            FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape(1)==11||peakshape(1)==12, % Fixed-width shapes only
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                if peakshape(1)==16||peakshape(1)==17, % Fixed-position shapes only
                    FitResults=[FitResults ; [round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]];
                else
                    FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]];
                end
            end
        end % if peakshape
    end % m==1
end % for m=1:NumPeaks
% Display Fit Results on upper graph
if plots,
    subplot(2,1,1);
    startx=min(xx)+xoffset+(max(xx)-min(xx))./20;
    dxx=(max(xx)-min(xx))./10;
    dyy=(max(yy)-min(yy))./10;
    starty=max(yy)-dyy;
    FigureSize=get(gcf,'Position');
    if peakshape(1)==9||peakshape(1)==10||peakshape(1)==19,  % Pulse and sigmoid shapes only
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
end % if plots

if NumArgOut==7,
    if plots,disp('Computing bootstrap sampling statistics.....'),end
    BootstrapResultsMatrix=zeros(5,100,NumPeaks);
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
        [FitResults,BootFitError]=fitpeaks(bx,by,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,FIXEDWIDTH);
        for peak=1:NumPeaks,
            BootstrapResultsMatrix(:,trial,peak)=FitResults(peak,:);
            BootstrapErrorMatrix(:,trial,peak)=BootFitError;
        end
    end
    if plots,toc;end
    for peak=1:NumPeaks,
        if plots,
            disp(' ')
            disp(['Peak #',num2str(peak) '         Position    Height       Width       Area']);
        end % if plots
        BootstrapMean=mean(real(BootstrapResultsMatrix(:,:,peak)'));
        BootstrapSTD=std(BootstrapResultsMatrix(:,:,peak)');
        BootstrapIQR=iqr(BootstrapResultsMatrix(:,:,peak)');
        PercentRSD=100.*BootstrapSTD./BootstrapMean;
        PercentIQR=100.*BootstrapIQR./BootstrapMean;
        BootstrapMean=BootstrapMean(2:5);
        BootstrapSTD=BootstrapSTD(2:5);
        BootstrapIQR=BootstrapIQR(2:5);
        PercentRSD=PercentRSD(2:5);
        PercentIQR=PercentIQR(2:5);
        if plots,
            disp(['Bootstrap Mean: ', num2str(BootstrapMean)])
            disp(['Bootstrap STD:  ', num2str(BootstrapSTD)])
            disp(['Bootstrap IQR:  ', num2str(BootstrapIQR)])
            disp(['Percent RSD:    ', num2str(PercentRSD)])
            disp(['Percent IQR:    ', num2str(PercentIQR)])
        end % if plots
        BootResults(peak,:)=[BootstrapMean BootstrapSTD PercentRSD BootstrapIQR PercentIQR];
    end % peak=1:NumPeaks,
end % if NumArgOut==6,
% ----------------------------------------------------------------------
function [FitResults,LowestError]=fitpeaks(xx,yy,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,fixedparameters)
% Based on peakfit Version 3: June, 2012. 
global PEAKHEIGHTS FIXEDWIDTH FIXEDPOSITIONS BIPOLAR
format short g
format compact
warning off all
FIXEDWIDTH=fixedparameters;
xoffset=0;
if start==0;start=calcstart(xx,NumPeaks,xoffset);end
PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
newstart=start;

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
        case 2
            TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart,options);
        case 3
            TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart,options);
        case 4
            TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart,options);
        case 5
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart,options);
        case 6
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart,options);
        case 7
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitlorentziancw(lambda,xx,yy)),cwnewstart,options);
        case 8
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart,options);
        case 9
            TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart,options);
        case 10
            TrialParameters=fminsearch(@(lambda)(fitupsigmoid(lambda,xx,yy)),newstar,optionst);
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
        case 14
            TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart,options);
        case 15
            TrialParameters=fminsearch(@(lambda)(fitBiLorentzian(lambda,xx,yy,extra)),newstart,options);
        case 16
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart,options);
        case 17
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart,options);
        case 18
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexplorentzian(lambda,zxx,zyy,-extra)),newstart,options);
        case 19
            TrialParameters=fminsearch(@(lambda)(alphafunction(lambda,xx,yy)),newstart,options);
        case 20
            TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart,options);
        case 21
            TrialParameters=fminsearch(@(lambda)(fittriangular(lambda,xx,yy)),newstart,options);
        case 22
            TrialParameters=fminsearch(@(lambda)(fitmultiple(lambda,xx,yy,shapesvector,extra)),newstart,options);
        case 23
            TrialParameters=fminsearch(@(lambda)(fitdownsigmoid(lambda,xx,yy)),newstar,optionst);
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
                A(m,:)=gaussian(xx,TrialParameters(m),FIXEDWIDTH);
            case 12
                A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDWIDTH);
            case 13
                A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 14
                A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 15
                A(m,:)=BiLorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
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
        end % switch
        %     for parameter=1:2:2*NumPeaks,
        %         newstart(parameter)=newstart(parameter)*(1+(rand-.5)/50);
        %         newstart(parameter+1)=newstart(parameter+1)*(1+(rand-.5)/10);
        %     end
    end % for
    
    % Multiplies each row by the corresponding amplitude and adds them up
    if AUTOZERO==3,
        baseline=PEAKHEIGHTS(1);
        Heights=PEAKHEIGHTS(2:1+NumPeaks);
        model=Heights'*A+baseline;
    else
        model=PEAKHEIGHTS'*A;
        Heights=PEAKHEIGHTS;
        baseline=0;
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

for m=1:NumPeaks,
    area(m)=trapz(xx+xoffset,height(m)*A(m,:)); % Compute the area of each component peak using trapezoidal method
end

for m=1:NumPeaks,
    if m==1,
        if peakshape(1)==6||peakshape(1)==7||peakshape(1)==8, % equal-width peak models
            FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape(1)==11||peakshape(1)==12,  % Fixed-width shapes only
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                FitResults=[[round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]];
            end
        end % if peakshape
    else
        if peakshape(1)==6||peakshape(1)==7||peakshape(1)==8, % equal-width peak models
            FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape(1)==11||peakshape(1)==12, % Fixed-width shapes only
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]];
            end
        end % if peakshape
    end % m==1
end % for m=1:NumPeaks
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
global PEAKHEIGHTS AUTOZERO BIPOLAR
numpeaks=round(length(lambda)/2);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(2*j-1),lambda(2*j))';
end 
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitewgaussian(lambda,t,y)
% Fitting function for a Gaussian band signal with equal peak widths.
global PEAKHEIGHTS AUTOZERO BIPOLAR
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),lambda(numpeaks+1))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFWGaussian(lambda,t,y)
%	Fitting function for a fixed width Gaussian
global PEAKHEIGHTS AUTOZERO FIXEDWIDTH BIPOLAR
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),FIXEDWIDTH)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFPGaussian(lambda,t,y)
%	Fitting function for fixed-position Gaussians
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,FIXEDPARAMETERS(j), lambda(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
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
global PEAKHEIGHTS AUTOZERO FIXEDWIDTH BIPOLAR
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),FIXEDWIDTH)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitewlorentzian(lambda,t,y)
% Fitting function for a Lorentzian band signal with equal peak widths.
global PEAKHEIGHTS AUTOZERO BIPOLAR
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),lambda(numpeaks+1))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
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
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = lorentzian(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
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
%	Fitting function for Logistic distribution, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  
%  T. C. O'Haver, May 2006
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = logistic(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = logistic(x,pos,wid)
% logistic distribution.  pos=position; wid=half-width (both scalar)
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
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = triangular(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = triangular(x,pos,wid)
%triangle function.  pos=position; wid=half-width (both scalar)
%triangle(x,pos,wid), where x may be scalar or vector,
%pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991
% Example
% x=[0:.1:10];plot(x,triangle(x,5.5,2.3),'.')
g=1-(1./wid) .*abs(x-pos);
for i=1:length(x),  
if g(i)<0,g(i)=0;end
end
% ----------------------------------------------------------------------
function err = fitpearson(lambda,t,y,shapeconstant)
%   Fitting functions for a Pearson 7 band signal.
% T. C. O'Haver (toh@umd.edu),   Version 1.3, October 23, 2006.
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = pearson(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
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
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = expgaussian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitexplorentzian(lambda,t,y,timeconstant)
%   Fitting functions for a exponentially-broadened lorentzian band signal.
%  T. C. O'Haver, 2013.
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = explorentzian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitexpewgaussian(lambda,t,y,timeconstant)
% Fitting function for exponentially-broadened Gaussian bands with equal peak widths.
global PEAKHEIGHTS AUTOZERO BIPOLAR
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = expgaussian(t,lambda(j),lambda(numpeaks+1),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
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
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = exppulse(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
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
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = alphafunction(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
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
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = downsigmoid(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitupsigmoid(tau,x,y)
% Fitting function for iterative fit to the sum of
% upwards moving sigmiods
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = upsigmoid(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g=downsigmoid(x,t1,t2)
% g=1./(1 + exp((t1 - x)./t2))'; % old version of sigmoid
g=.5-.5*erf(real((x-t1)/sqrt(2*t2))); % down step sigmoid
% ----------------------------------------------------------------------
function g=upsigmoid(x,t1,t2)
% g=1./(1 + exp((t1 - x)./t2))'; % old version of sigmoid
g=1/2 + 1/2* erf(real((x-t1)/sqrt(2*t2))); % up step sigmoid
% ----------------------------------------------------------------------
function err = fitGL(lambda,t,y,shapeconstant)
%   Fitting functions for Gaussian/Lorentzian blend.
% T. C. O'Haver (toh@umd.edu), 2012.
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = GL(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
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
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = voigt(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
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
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = BiGaussian(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
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
function err = fitBiLorentzian(lambda,t,y,shapeconstant)
%   Fitting functions for BiGaussian.
% T. C. O'Haver (toh@umd.edu),  2012.
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = BiLorentzian(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = BiLorentzian(x,pos,wid,m)
% BiLorentzian (different widths on leading edge and trailing edge).
% pos=position; wid=width 
% m determines shape; symmetrical if m=50.
%  T. C. O'Haver, 2012
lx=length(x);
hx=val2ind(x,pos);
g(1:hx)=lorentzian(x(1:hx),pos,wid*(m/100));
g(hx+1:lx)=lorentzian(x(hx+1:lx),pos,wid*(1-m/100));
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
global PEAKHEIGHTS AUTOZERO BIPOLAR
numpeaks=round(length(lambda)/2);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = peakfunction(shapesvector(j),t,lambda(2*j-1),lambda(2*j),m(j))';
end 
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function p=peakfunction(shape,x,pos,wid,m)
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
        p=expgaussian(x,pos,wid,m)';
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
        p=BiLorentzian(x,pos,wid,m);
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
        p=downpsigmoid(x,pos,wid);
    otherwise
end % switch
% ----------------------------------------------------------------------
function ShapeString=SelectShapeString(Shape)
global shapesvector
switch Shape
    case 1
        ShapeString='Gaussian';
    case 2
        ShapeString='Lorentzian';
    case 3
        ShapeString='Logistic distribution';
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
        ShapeString='bifurcated Lorentzian';
    case 16
        ShapeString='Fixed-position Gaussians';
    case 17
        ShapeString='Fixed-position Lorentzians';
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
    otherwise
        ShapeString='';
end % switch Shape
% ----------------------------------------------------------------------
