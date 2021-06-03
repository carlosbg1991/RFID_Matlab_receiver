function pY=isignal(DataMatrix,xcenter,xrange,sm,sw,em,dm,rm,s1,s2,sr,mw)
% Y=isignal(DataMatrix,xcenter,xrange,SmoothMode,SmoothWidth,ends,...
% DerivativeMode,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth)
% Version 2.5. See http://terpconnect.umd.edu/~toh/spectrum/iSignal.html
% An interactive function that performs smoothing, differentiation, and
% peak sharpening of a time-series signal in the form of a 2-column 
% matrix with the independent variable (x-values) in the first
% column and dependent variable (y values) in the second column, or as
% separate x and y vectors. Returns the processed independent axis (Y)
% vector as the output argument. The lower half of the figure window shows
% a plot of the entire signal, and the upper half shows a selected portion
% controlled by the pan and zoom keystrokes or by optional input  
% arguments 'xcenter' and 'xrange', respectively. Other keystrokes
% also allow you to control the smooth type, width, and ends
% treatment, the derivative order (0th through 5th), and peak 
% sharpening. (Alternatively, the initial values of these parameters 
% can be passed to the function via the optional input arguments 
% sm,sw,em,dm,rm,s1,s2). 
% Version 2.2: Adds slew rate limit (~ key) and spike filter (M key).
% Version 2.3: Measures peak area two ways. Fixes several bugs related to input arguments.
% Version 2.4: Computes signal-to-noise ratio (SNR) in Peak mode.
% Version 2.5: Bug fixes
% By T. C. O'Haver (toh@umd.edu); Savitzky-Golay smooth code by Diederick.  
%
% The S key (or optional argument "sm") determines the smooth mode:
%      If sm=0, the signal is not smoothed.
%      If sm=1, rectangular (sliding-average or boxcar) 
%      If sm=2, triangular (2 passes of sliding-average)
%      If sm=3, pseudo-Gaussian (3 passes of sliding-average)
%      If sm=4, Savitzky-Golay smooth 
% The A and Z keys (or optional argument sw) control the smooth width.
% The Z key (or argument "em") controls how the "ends" of the signal 
%    (the first w/2 points and the last w/2 points) are handled.
%      If ends=0, the ends are zeroed
%      If ends=1, the ends are smoothed with progressively 
%      smaller smooths the closer to the end.
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html
%
% The D key (or optional input argument "dm") determines the derivative 
%    order (O, 1, 2, 3, 4, 5, and back to 0). See
%    http://terpconnect.umd.edu/~toh/spectrum/Differentiation.html
%
% The E key (or optional argument "rm") turns off and on peak 
% sharpening (resolution enhancement). The sharpening strength is
% controled by the F and V keys (optional argument "s1") and B and G 
% keys (optional argument "s2"). The optimum values depend on the 
% peak shape and width; For details, see
% http://terpconnect.umd.edu/~toh/spectrum/InteractiveResEnhance.htm). 
%
% The P key toggles off and on the peak measure mode, which measures and
% displays the peak position, height, width, and area of the one peak
% at a time if it is centered and zoomed in; a red "cap" on the peak
% indicates that portion of the signal that is taken for the measurement.
% Press 'R' key to print out the peak measures in the command window.
% The L key toggles off and on the Overlay mode, which overlays the
% selected portion in the upper plot with the original signal as a 
% dotted line, for comparison. The H key switches between linear and log
% y-axis on the lower plot. The 0 (zero) key set minimun signal to zero.
% The ; key (semicolon) sets the entire selected region to zero (use to
% remove stray data points). The Tab key resets smooth, derivative, and 
% sharpen effects to zero.  The O (letter O) key saves the X,Y processed 
% signal as a "mat" file, in a location and with a file name that you 
% specify. The C key condenses the signal by the specified factor N,
% replacing each group of N points with their average; the I key replaces 
% the signal with a linearily interploated version containing N data points.
% The M key implements a median filter for removing spikes. The ~ key
% limits the maximum rate of change of the signal.
%
% Press K to see all keyboard commands.
%
% EXAMPLE 1: Data in two columns of a matrix: [x y].
%             >> load data.mat
%             >> isignal(DataMatrix);
% 
% EXAMPLE 2: Data in separate x,y vectors or single y vector
%             >> isignal(x,y);  or
%             >> isignal(y);  
%
% EXAMPLE 3: As above, but specifies initial values of pan (xcenter) and 
%            zoom (xrange) in the last two input arguments. 
%             >> isignal(DataMatrix,180,40); or
%             >> isignal(x,y,180,40);
%
% EXAMPLE 4: As above, but additionally specifies initial values of 
%            SmoothMode, SmoothWidth, ends, and DerivativeMode. 
%             >> isignal(DataMatrix,180,40,2,9,0,1);
% 
% EXAMPLE 5: As above, but additionally specifies initial values of the  
%            peak sharpening parameters Sharpen, Sharp1, and Sharp2.
%             >> isignal(DataMatrix,180,40,4,19,0,0,1,51,6000);
%                (Press 'E' key to toggle sharpening ON/OFF)
%
% EXAMPLE 6:   >> x=[0:.005:2];y=humps(x);Data=[x;y];
%             4th derivative of the peak at x=0.9:
%              >> isignal(Data,0.9,0.5,1,3,1,4);
%             Peak sharpening applied to the peak at x=0.3:
%              >> isignal(Data,0.3,0.5,4,3,1,0,1,220,5400);
%                 (Press 'E' key to toggle sharpening ON/OFF)
%
% EXAMPLE 7: Measurement of peak area.  This example generates four 
% Gaussian peaks, all with the exact same peak height (1.00) and area 
% (1.77). The first peak (at x=4) is isolated, the second peak (x=9) 
% is slightly overlapped with the third one, and the last two peaks 
% (at x= 13 and 15) are strongly overlapped.  To measure the area under 
% a peak using the perpendicular drop method, position the dotted red
% marker lines at the minimum between the overlapped peaks.  
% 
% >>  x=[0:.01:20];y=exp(-(x-4).^2)+exp(-(x-9).^2)+exp(-(x-13).^2)+exp(-(x-15).^2);isignal(x,y);
%
% EXAMPLE 8: Single peak with random spikes. Compare smoothing vs spike
%            filter (M key) and slew rate limit (~ key) to remove spikes.
% >> x=-5:.01:5;
% >> y=exp(-(x).^2);for n=1:1000,if randn()>2,y(n)=rand()+y(n);,end,end;
% >> isignal(x,y);
%
% KEYBOARD CONTROLS:
%  Pan signal left and right...Coarse pan: < and >
%                              Fine pan: left and right cursor arrows
%                              Nudge: [ and ]
%  Zoom in and out.............Coarse zoom: / and "  
%                              Fine zoom: up and down cursor arrows
%  Resets pan and zoom.........ESC
%  Adjust smooth width.........A,Z  (A=>more, Z=>less) 
%  Adjust smooth type..........S (No, Rectanglular, Triangle, Gaussian, Savitzky-Golay)
%  Toggle smooth ends..........X (0=ends zeroed  1=ends smoothed (slower)
%  Adjust derivative order.....D (0th to 5th order)
%  Toggle peak sharpening......E (0=OFF 1=ON)
%  Sharpening for Gaussian.....Y  Set sharpen settings for Gaussian
%  Sharpening for Lorentzian...U  Set sharpen settings for Lorentzian
%  Adjust sharp1...............F,V  F=>sharper, V=>less sharpening
%  Adjust sharp2...............G,B  G=>sharper, B=>less sharpening
%  Slew rate limit (0=OFF).....`  Largest allowed change between points
%  Spike filter width (O=OFF)..m  Spike filter eliminates sharp spikes
%  Toggle peak parabola........P  fits parabola to center, labels vertex
%  Print peak report...........R  prints position, height, width, area
%  Toggle overlay mode.........L  Overlays original signal as dotted line
%  Toggle log y mode...........H  semilog plot in lower window
%  Toggle Autozero OFF/ON......T  Subtracts background from upper window segment
%  Restores original signal....Tab key resets to original signal and modes
%  Baseline subtraction........Backspace, then click baseline at multiple points
%  Restore background..........\  to cancel previous background subtraction
%  Invert signal...............-  Invert (negate) the signal (flip + and -)
%  Remove offset...............0  (zero) set minimun signal to zero 
%  Sets region to zero.........;  sets selected region to zero.
%  Condense signal.............C  Condense oversampled signal by factor N
%  Interpolate signal..........i  Interpolate (resample) to N points
%  Print keyboard commands.....K  prints this list
%  Print signal report.........Q  prints signal info and current settings
%  Print isignal arguments.....W  prints isignal (current arguments)
%  Save output to disk.........O as .mat file with processed signal matrix

global X Y xo dx DerivativeMode Sharpen Sharp1 Sharp2 SmoothWidth SlewRate MedianWidth
global SmoothType ends SmoothMode SavedSignal SavedXvalues PeakLabels Report autozero
format short g
format compact
warning off all
switch nargin % Process arguments
    % 'nargin' is the number of arguments
    case 1  % One argument only
        % Might be isignal(DataMatrix) ot isignal(Y-vector)
        % If data is in the wrong transposition, fix it.
        datasize=size(DataMatrix);
        if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
        datasize=size(DataMatrix);
        if datasize(2)==1, %  Must be isignal(Y-vector)
            X=1:length(DataMatrix); % Create an independent variable vector
            Y=DataMatrix;
        else
            % Must be isignal(DataMatrix)
            X=DataMatrix(:,1); % Split matrix argument
            Y=DataMatrix(:,2);
        end
        SmoothMode=0; % Initial SmoothMode = Rect
        SmoothWidth=3; % Initial smooth width
        SmoothType='Rect.';  % Label for initial SmoothMode
        DerivativeMode=0; % Derivative mode off initially
        ends=0;  % Initial smooth ends setting zero (no tapering)
        Sharpen=0;  % Initially no sharpening
        Sharp1=10; % Initial factor1 for resolution enhancement
        Sharp2=1000; % Initial factor2 for resolution enhancement
        SlewRate=0;
        MedianWidth=0;
        xo=length(Y)/2; % Initial Pan setting
        dx=length(Y)/4; % Initial Zoom setting
    case 2
        % Two arguments, might be separate x and y data vectors,
        % or one data matrix and a peak density estimate.
        if isscalar(xcenter) % if second argument is scalar
            % Must be isignal(DataMatrix,xcenter)
            % If DataMatrix is in the wrong transposition, fix it.
            datasize=size(DataMatrix);
            if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
            X=DataMatrix(:,1); % Split matrix argument
            Y=DataMatrix(:,2);
            xo=val2ind(X,xcenter);
        else % if second argument is not scalar
            % Must be isignal(x,y)
            xdatasize=size(DataMatrix);
            if xdatasize(1)<xdatasize(2),DataMatrix=DataMatrix';end
            X=DataMatrix;  % First argument is X
            xdatasize=size(xcenter);
            if xdatasize(1)<xdatasize(2),xcenter=xcenter';end
            Y=xcenter; % Second argument is Y
            xo=length(Y)/2; %  % Default initial zoom setting
        end  % if isscalar
        SmoothMode=0; % Initial SmoothMode = No
        SmoothWidth=3; % Initial smooth width
        SmoothType='Rect.';  % Label for initial SmoothMode
        DerivativeMode=0; % Derivative mode off initially
        ends=0;  % Initial smooth ends setting zero (no tapering)
        Sharpen=0;  % Initially no sharpening
        Sharp1=10; % factor1 for resolution enhancement
        Sharp2=1000; % factor2 for resolution enhancement
         SlewRate=0;
         MedianWidth=0;
        dx=length(Y)/4; %  % Default initial zoom setting
    case 3
        % Might be isignal(DataMatrix,xcenter,xrange) or isignal(x,y,xcenter)
        if isscalar(xcenter) % if second argument is scalar
            % Must be isignal(DataMatrix,xcenter,xrange)
            % If DataMatrix is in the wrong transposition, fix it.
            datasize=size(DataMatrix);
            if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
            X=DataMatrix(:,1); % Split matrix argument
            Y=DataMatrix(:,2);
            [xo,dx]=panandzoom(X,xcenter,xrange);
        else % if second argument is not isscalar
            % Must be isignal(x,y,xcenter)
            xdatasize=size(DataMatrix);
            if xdatasize(1)<xdatasize(2),DataMatrix=DataMatrix';end
            X=DataMatrix;  % First argument is X
            xdatasize=size(xcenter);
            if xdatasize(1)<xdatasize(2),xcenter=xcenter';end
            Y=xcenter; % Second argument is Y
            xo=xrange; % third argument is xcenter
            dx=length(Y)/4; % Default initial zoom setting
        end  % if isscalar
        SmoothMode=0; % Initial SmoothMode = No
        SmoothWidth=3; % Initial smooth width
        SmoothType='Rect.';  % Label for initial SmoothMode
        DerivativeMode=0; % Derivative mode off initially
        ends=0;  % Initial smooth ends setting zero (no tapering)
        Sharpen=0;  % Initially no sharpening
        Sharp1=10; % factor1 for resolution enhancement
        Sharp2=1000; % factor2 for resolution enhancement
         SlewRate=0;
         MedianWidth=0;
    case 4   % Must be isignal(x,y,xcenter,xrange)
        xdatasize=size(DataMatrix);
        if xdatasize(1)<xdatasize(2),DataMatrix=DataMatrix';end
        X=DataMatrix;  % First argument is X
        xdatasize=size(xcenter);
        if xdatasize(1)<xdatasize(2),xcenter=xcenter';end
        Y=xcenter; % Second argument is Y
        SmoothMode=0; % Initial SmoothMode = No
        SmoothWidth=3; % Initial smooth width
        SmoothType='Rect.';  % Label for initial SmoothMode
        DerivativeMode=0; % Derivative mode off initially
        ends=0;  % Initial smooth ends setting zero (no tapering)
        Sharpen=0;  % Initially no sharpening
        Sharp1=10; % factor1 for resolution enhancement
        Sharp2=1000; % factor2 for resolution enhancement
         SlewRate=0;
         MedianWidth=0;
        [xo,dx]=panandzoom(X,xrange,sm);
    case 7
        % One data matrix, all smoothing and derivative parameters specified
        % in arguments, default values for resolution enhancement.
        % If DataMatrix is in the wrong transposition, fix it.
        datasize=size(DataMatrix);
        if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
        X=DataMatrix(:,1); % Split matrix argument
        Y=DataMatrix(:,2);
        [xo,dx]=panandzoom(X,xcenter,xrange);
        SmoothMode=sm; % SmoothMode (O, 1, 2, 3 or 4)
        SmoothWidth=sw; % Smooth width
        DerivativeMode=dm; % Derivative mode (0, 1, 2, 3, 4)
        ends=em;  % Smooth ends setting (0 or 1)
        Sharpen=0;  % Initially no sharpening
        Sharp1=10; % factor1 for resolution enhancement
        Sharp2=1000; % factor2 for resolution enhancement
         SlewRate=0;
         MedianWidth=0;
    case 10
        % One data matrix, all signal processing parameters specified
        % in arguments, including  resolution enhancement.
        % If DataMatrix is in the wrong transposition, fix it.
        datasize=size(DataMatrix);
        if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
        X=DataMatrix(:,1); % Split matrix argument
        Y=DataMatrix(:,2);
        [xo,dx]=panandzoom(X,xcenter,xrange);
        SmoothMode=sm; % SmoothMode sm
        SmoothWidth=sw; % Smooth width
        DerivativeMode=dm; % Derivative mode (0, 1, 2, 3, 4)
        ends=em;  % Smooth ends setting (0 or 1)
        Sharpen=rm;  % Sharpen mode
        Sharp1=s1; % factor1 for resolution enhancement
        Sharp2=s2; % factor2 for resolution enhancement
         SlewRate=0;
         MedianWidth=0;
    case 11
        % One data matrix, all signal processing parameters specified
        % in arguments, including  resolution enhancement.
        % If DataMatrix is in the wrong transposition, fix it.
        datasize=size(DataMatrix);
        if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
        X=DataMatrix(:,1); % Split matrix argument
        Y=DataMatrix(:,2);
        [xo,dx]=panandzoom(X,xcenter,xrange);
        SmoothMode=sm; % SmoothMode sm
        SmoothWidth=sw; % Smooth width
        DerivativeMode=dm; % Derivative mode (0, 1, 2, 3, 4)
        ends=em;  % Smooth ends setting (0 or 1)
        Sharpen=rm;  % Sharpen mode
        Sharp1=s1; % factor1 for resolution enhancement
        Sharp2=s2; % factor2 for resolution enhancement
         SlewRate=sr;
         MedianWidth=0;
    case 12
        % One data matrix, all signal processing parameters specified
        % in arguments, including  resolution enhancement.
        % If DataMatrix is in the wrong transposition, fix it.
        datasize=size(DataMatrix);
        if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
        X=DataMatrix(:,1); % Split matrix argument
        Y=DataMatrix(:,2);
        [xo,dx]=panandzoom(X,xcenter,xrange);
        SmoothMode=sm; % SmoothMode sm
        SmoothWidth=sw; % Smooth width
        DerivativeMode=dm; % Derivative mode (0, 1, 2, 3, 4)
        ends=em;  % Smooth ends setting (0 or 1)
        Sharpen=rm;  % Sharpen mode
        Sharp1=s1; % factor1 for resolution enhancement
        Sharp2=s2; % factor2 for resolution enhancement
        SlewRate=sr;
        MedianWidth=mw;
    otherwise
        disp('Invalid number of arguments')
        disp('Expected forms are:')
        disp('isignal(y);  % Data in single y vector')
        disp('isignal(x,y);  % Data in separate x and y vectors')
        disp('isignal(DataMatrix); % Data in two columns of DataMatrix')
        disp('isignal(x,y,xcenter,xrange); ')
        disp('isignal(DataMatrix,xcenter,xrange;) ')
        disp('isignal(DataMatrix,xcenter,xrange,SmoothMode,SmoothWidth,ends,DerivativeMode); ' )
        disp('isignal(DataMatrix,xcenter,xrange,SmoothMode,SmoothWidth,ends,DerivativeMode,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth); ')
        beep
        return
end % switch nargin
% Define smooth type string for xlabel
switch SmoothMode
     case 0
          SmoothType='No';
     case 1
          SmoothType='Rect.';
     case 2
          SmoothType='Tri.';
     case 3
          SmoothType='Gauss';
     case 4
           SmoothType='Savitzky-Golay';
end
PeakLabels=0;  % Start with peak label turned off
% Save original signal in SavedSignal for undo function 
SavedSignal=Y;
SavedXvalues=X;
Overlay=0;  % Start with overlay turned off
Report=0;
autozero=0;
logymode=0; % Start in linear y mode.

% Plot the signal
pY=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
Y=pY;
RedrawSignal(X,Y,xo,dx);

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
% add that to the SWITCH statement for your own custom functions.
global X Y xx yy xo dx SmoothMode SmoothWidth DerivativeMode SlewRate
global Sharpen Sharp1 Sharp2 SavedSignal SavedXvalues SavedBackground
global SmoothType ends PeakLabels Overlay Report GaussEstimate MedianWidth
global LorentzEstimate logymode autozero
key=get(gcf,'CurrentCharacter');
if isscalar(key),
    ly=length(Y);
    switch double(key),
        case 28
            % Pans down when right arrow pressed.
            xo=xo-dx/20;
            if xo<1,xo=1;end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 29
            % Pans up when left arrow pressed.
            xo=xo+dx/20;
            if xo>ly,xo=ly;end
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
        case 44
            % Pans down when > key pressed.
            xo=xo-dx;
            if xo<1,xo=1;end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 46
            % Pans up when < key pressed.
            xo=xo+dx;
            if xo>ly,xo=ly;end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 31
            % Zooms up when up arrow pressed.
            dx=dx+dx/10;
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 30
            % Zooms down when down arrow pressed.
            dx=dx-dx/10;
            if dx<2,dx=2;end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 47
            % Zooms x 2 up when / pressed.
            dx=dx*2;
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 39
            % Zooms x 1/2 down when ' pressed.
            dx=round(dx/2);
            if dx<2,dx=2;end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 8
            % When 'Backspace' key is pressed, user clicks the graph
            % along the presumed background points, then the program
            % subtracts the interploated background between the points.
            SavedBackground=Y;
            disp('Multi-point baseline subtraction')
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
            tempY=Y;
            for k=1:length(bX)-1,
                fp=val2ind(X,bX(k)); % First point in segment
                lp=val2ind(X,bX(k+1));  % Last point in segment
                % Subtract piecewise linear background from Y
                tempY(fp:lp)=Y(fp:lp)-((bY(k+1)-bY(k))/(bX(k+1)-bX(k))*(X(fp:lp)-bX(k))+bY(k));
            end
            Y=tempY;
            SavedSignal=Y;
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 92
            % When '\' key is pressed, restoreed original signal
            SavedSignal=SavedBackground;
            X=SavedXvalues;
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 122
            % When 'z' key is pressed, DEcreases "SmoothWidth" by 1 or 10%
            if SmoothMode==0,
                SmoothMode=1;
                SmoothType='Rect.';
            end
            if SmoothWidth>20,
                SmoothWidth=round(SmoothWidth-.1.*SmoothWidth);
                SmoothWidth=2*round(SmoothWidth/2)-1;
            else
                SmoothWidth=SmoothWidth-2;
            end
            if SmoothWidth<1, SmoothWidth=1;end
            if SmoothMode==0,
                Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            else
                Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 97
            % When 'a' key is pressed, INcreases "SmoothWidth" by 1 or 10%
            if SmoothMode==0,
                SmoothMode=1;
                SmoothType='Rect.';
            end
            if SmoothWidth>20,
                SmoothWidth=round(SmoothWidth+.1.*SmoothWidth);
                SmoothWidth=2*round(SmoothWidth/2)+1;
            else
                SmoothWidth=SmoothWidth+2;
            end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 115 % When 's' key is pressed, steps through SmoothModes
            SmoothMode=SmoothMode+1;
            if SmoothMode==5,SmoothMode=0; end
            switch SmoothMode
                case 0
                    SmoothType='No';
                case 1
                    SmoothType='Rect.';
                case 2
                    SmoothType='Tri.';
                case 3
                    SmoothType='Gauss';
                case 4
                    SmoothType='Savitzky-Golay';
            end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 120 % When 'x' key is pressed, toggles between ends 0 and 1
            if ends==0,
                ends=1;
            else
                ends=0;
            end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 100
            % When 'd' key is pressed, cycles through DerivativeModes 0,1,2,3,4,5->0
            % if length(Y)>10000,disp('Warning: Derivatives can be slow for for signal lengths above 10,000 points'),end
            DerivativeMode=DerivativeMode+1;
            if DerivativeMode==6,DerivativeMode=0; end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 101 % When 'e' key is pressed, toggles between Sharpen 0 and 1
            % if length(Y)>10000,disp('Warning: Sharpening can be slow for for signal lengths above 10,000 points'),end
            if Sharpen==0,
                Sharpen=1;
                SmoothMode=4;
                if SmoothWidth<3;SmoothWidth=3;end
                SmoothType='Savitzky-Golay';
            else
                Sharpen=0;
            end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 121  % When 'y' key is pressed, sets Sharp1 and 2 for Gaussian
            GaussEstimate=1;
            % if length(Y)>10000,disp('Warning: Sharpening can be slow for
            % for signal lengths above 10,000 points'),end
            PeakLabels=1;
            SmoothMode=4;
            if SmoothWidth<3;SmoothWidth=3;end
            SmoothType='Savitzky-Golay';
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
            Sharpen=1;
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 117  % When 'u' key is pressed, sets Sharp1 and 2 for Lorentzian
            LorentzEstimate=1;
            % if length(Y)>10000,disp('Warning: Sharpening can be slow for for signal lengths above 10,000 points'),end
            SmoothMode=4;
            if SmoothWidth<3;SmoothWidth=3;end
            SmoothType='Savitzky-Golay';
            PeakLabels=1;
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
            Sharpen=1;
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 102
            % When 'f' key is pressed, increases Sharp1
            if Sharpen==0,Sharpen=1;end
            Sharp1=Sharp1+.1*Sharp1;
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 118
            % When 'v' key is pressed, decreases Sharp1
            if Sharpen==0,Sharpen=1;end
            Sharp1=Sharp1-.1*Sharp1;
            if Sharp1<0, Sharp1=0;end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 103
            % When 'g' key is pressed, increases Sharp2
            if Sharpen==0,Sharpen=1;end
            Sharp2=Sharp2+.1*Sharp2;
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 98
            % When 'b' key is pressed, decreases Sharp2
            if Sharpen==0,Sharpen=1;end
            Sharp2=Sharp2-.1*Sharp2;
            if Sharp2<0, Sharp2=0;end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 112
            % When 'p' key is pressed, toggles on/off peak labels in upper panel
            if PeakLabels==0,
                PeakLabels=1;
            else
                PeakLabels=0;
            end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 108 % When 'L' key is pressed, toggles between ends 0 and 1
            if Overlay==0,
                Overlay=1;
            else
                Overlay=0;
            end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 104 % When 'H' key is pressed, toggles between normal and logy plot
            if logymode==0,
                logymode=1;
            else
                logymode=0;
            end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 114 % When 'R' key is pressed, toggles between Report 0 and 1
            if Report==0,
                Report=1;
                if autozero,disp('Autozero ON'),end
                disp('Position     Height       Width     Gaussian Area   Total Area     SNR')
            else
                Report=0;
            end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 116 % When 'T' key is pressed, toggles between normal and autozero plot
            if autozero==0,
                autozero=1;
            else
                autozero=0;
            end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 9 % When 'Tab' key is pressed, resets to original signal and modes
            Y=SavedSignal;
            X=SavedXvalues;
            SmoothMode=0;
            SmoothWidth=1;
            SmoothType='No';
            DerivativeMode=0;
            Sharpen=0;
            SlewRate=0;
            MedianWidth=0;
            Y=ProcessSignal(X,Y,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 27 % When 'ESC' key is pressed, resets pan and zoom
            xo=length(Y)/2; % Initial Pan setting
            dx=length(Y)/4; % Initial Zoom setting
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 59
            % When ';' key is pressed, replaces selected segment with zeros
            doit=input('Replace selected region with zeros?','s');
            if doit=='y',
                startpoint=round(xo-dx/2);
                if startpoint<1;startpoint=1;end
                endpoint=round(xo+dx/2);
                if endpoint>length(Y);endpoint=length(Y);end
                Y(startpoint:endpoint)=0;
            end
             [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 48
            % When '0' key is pressed, minimum value of Y is set to zero
            Y=Y-min(Y);
            SavedSignal=Y;
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
            disp('Signal re-zeroed')
        case 45
            % When '-' (dash or minus) key is pressed, invert the signal
            SavedSignal=-SavedSignal;
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 96
            % When '~' (tilde) key is pressed, enforce maximum slew rate
            SavedSignal=Y;
            disp(['Current slew rate limit =' num2str(SlewRate)])
            SlewRate=input('Enter desired slew rate:');
            if SlewRate=='',SlewRate=0;end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 109
            % When 'm' key is pressed, performs median filter
            SavedSignal=Y;
            disp(['Current spike width =' num2str(MedianWidth)])
            MedianWidth=input('Enter spike width (1,2,3,...):');
            if MedianWidth=='',MedianWidth=0;end
            MedianWidth=round(MedianWidth);
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 127
            % When 'Delete' key is pressed, sets the single point under the
            % green cursor to zero
            Y(round(xo))=0;
            SavedSignal=Y;
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
%       case 95
%             % When '_' key (Shift '-') is pressed, replaces selected region
%             startpoint=round(xo-dx/2);
%             if startpoint<1;startpoint=1;end
%             endpoint=round(xo+dx/2)-1;
%             if endpoint>length(Y);endpoint=length(Y);end
%             lxx=length(xx);
%             bkgsize=2;
%             X1=xx(1:round(lxx/bkgsize));
%             X2=xx((lxx-round(lxx/bkgsize)):lxx);
%             MX=[X1;X2];
%             Y1=yy(1:round(length(xx)/bkgsize));
%             Y2=yy((lxx-round(lxx/bkgsize)):lxx);
%             MY=[Y1;Y2];
%             bkgcoef=polyfit(MX,MY,1);  % Fit straight line to sub-group of points
%             bkg=polyval(bkgcoef,xx);
%             Y(startpoint:endpoint)=bkg;
%             Y=ProcessSignal(X,Y,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
%             [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 73
            % When 'I' key (upper-case i or Shift-i) is pressed, integrates the signal
            sum=0;
            for n=1:length(X),
                sum=sum+Y(n);
                Y(n)=sum;
            end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 105
            % When 'i' key (lower-case i) is pressed, interpolates the signal
            % to find XI,YI, the values of the underlying function Y at the points
            % linearly interpolated between the points of X.
            disp(['X,Y size before interpolation = ' num2str(size(X)) ' , '  num2str(size(Y)) ] )
            InterPoints=input('Number of points in interpolated signal: ');
            if InterPoints>1,
                Xi=linspace(min(X),max(X),InterPoints);
                Y=interp1(X,Y,Xi)';
                X=Xi';
                disp(['X,Y size after interpolation = ' num2str(size(X)) ' , '  num2str(size(Y)) ] )
                xo=length(Y)/2; % Initial Pan setting
                dx=length(Y)/4; % Initial Zoom setting
                SavedSignal=Y;
                SavedXvalues=X;
                pY=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
                Y=pY;
                RedrawSignal(X,Y,xo,dx);
            end
        case 99
            % When C key is pressed, condenses signal by specified factor
            CondenseFactor=input('Condense oversampled signal by factor of (e.g. 2, 3, 4...): ');
            if CondenseFactor>1,
                disp([ 'X,Y size before condensation = ' num2str(size(X)) ' , '  num2str(size(Y)) ] )
                X=condense(X,CondenseFactor)';
                Y=condense(Y,CondenseFactor)';
                xo=length(Y)/2; % Initial Pan setting
                dx=length(Y)/4; % Initial Zoom setting
                SavedSignal=Y;
                SavedXvalues=X;
                disp([ 'X,Y size after condensation = ' num2str(size(X)) ' , '  num2str(size(Y)) ] )
                pY=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
                Y=pY;
                RedrawSignal(X,Y,xo,dx);
            end
        case 113
            % When 'Q' is pressed, prints a report listing signal
            % characteristics and current settings.
            disp('--------------------------------------------------------')
            % SizeX=size(X)
            % SizeY=size(Y)
            disp(['X,Y size = ' num2str(size(X)) ' , '  num2str(size(Y)) ] )
            disp([ num2str(length(Y)) ' total points from x= ' num2str(X(1)) ' to '  num2str(X(length(X))) ] )
            disp(['Interval between points = ' num2str(X(2)-X(1)) ' to ' num2str(X(length(X))-X(length(X)-1)) ] )
            disp(sprintf('min/max Y = %0.3g / %0.4g', min(Y), max(Y)))
            if autozero,
                disp('Autozero ON'),
            else
                disp('Autozero OFF'),
            end
            if SlewRate,
                disp(['Maximum slew rate = ' num2str(SlewRate) ] ),
            end
            if MedianWidth,
                disp(['spike filter width = ' num2str(MedianWidth) ] ),
            end
            disp(['Smooth: = ' num2str(SmoothWidth) ' point ' SmoothType ', Ends = ' num2str(ends) ] )
            if DerivativeMode,
                disp(['Derivative order = ' num2str(DerivativeMode) ] ),
            end
            if Sharpen
                disp(['Sharpen factor 1 = ' num2str(Sharp1) '  Sharpen factor 2 = ' num2str(Sharp2) ] )
            end
            disp([ 'Selected range: ' num2str(length(xx)) ' points from x=' num2str(min(xx)) ' to ' num2str(max(xx)) ])
            disp(sprintf('  Peak-to-peak Y: %0.4g  \r  Standard deviation: %0.3g ', max(yy)-min(yy), std(yy)))
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
            disp(sprintf('  Area: %0.4g', trapz(xx,yy)))
        case 107
            % When 'k' key is pressed, prints out table of keyboard commands
            disp('KEYBOARD CONTROLS, version 2.5:')
            disp(' Pan signal left and right...Coarse pan: < and >')
            disp('                             Fine pan: left and right cursor arrows')
            disp('                             Nudge: [ and ] ')
            disp(' Zoom in and out.............Coarse zoom: / and "  ')
            disp('                             Fine zoom: up and down cursor arrows')
            disp(' Resets pan and zoom.........ESC')
            disp(' Adjust smooth width.........A,Z (A=>more, Z=>less) ')
            disp(' Cycle smooth types..........S (No, Rectanglular, Triangle, Gaussian, Savitzky-Golay)')
            disp(' Toggle smooth ends..........X (0=ends zeroed  1=ends smoothed (slower)')
            disp(' Cycle derivative orders.....D (Steps through 0th to 5th order)')
            disp(' Toggle peak sharpening......E (0=OFF 1=ON)')
            disp(' Sharpening for Gaussian.....Y  Set sharpen settings for Gaussian')
            disp(' Sharpening for Lorentzian...U  Set sharpen settings for Lorentzian')
            disp(' Adjust sharp1...............F,V  F=>sharper, V=>less sharpening')
            disp(' Adjust sharp2   ............G,B  G=>sharper, B=>less sharpening')
            disp(' Slew rate limit (0=OFF).....~  Largest allowed y change between points')
            disp(' Spike filter width (O=OFF)..m  spike filter eliminates sharp spikes')
            disp(' Toggle peak parabola........P  fits parabola to center, labels vertex')
            disp(' Print peak report...........R  prints position, height, width, area')
            disp(' Toggle log y mode...........H  semilog plot in lower window')
            disp(' Toggle Autozero OFF/ON......T  Subtracts background from upper window segment')
            disp(' Restores original signal....Tab key resets to original signal and modes')
            disp(' Toggle overlay mode.........L  Overlays original signal as dotted line')
            disp(' Baseline subtraction........Backspace, then click baseline at multiple points')
            disp(' Restore background..........\  to cancel previous background subtraction')
            disp(' Invert signal...............-  Invert (negate) the signal (flip + and -)')
            disp(' Remove offset...............0  (zero) set minimun signal to zero ')
            disp(' Sets region to zero.........;  sets selected region to zero')
            disp(' Condense signal.............C  Condense oversampled signal by factor of N')
            disp(' Interpolate signal..........i  Interpolate (resample) to N points')
            disp(' Print report................Q  prints signal info and current settings')
            disp(' Print keyboard commands.....K  prints this list')
            disp(' Print isignal arguments.....W  prints isignal function with all current arguments')
            disp(' Save output to disk.........O  Save .mat file with processed signal matrix')
        case 119
            % When 'W' is pressed, prints 'isignal(current arguments)'
            firstpoint=xo+dx/2;
            if firstpoint>length(X),firstpoint=length(X);end
            lastpoint=xo-dx/2;
            if lastpoint<1,lastpoint=1;end
            disp(['isignal(DataMatrix,'  num2str(X(round(xo))) ',' num2str(X(round(firstpoint))-X(round(lastpoint)))  ',' num2str(SmoothMode)  ',' num2str(SmoothWidth) ',' num2str(ends) ',' num2str(DerivativeMode)  ',' num2str(Sharpen)  ',' num2str(Sharp1)  ',' num2str(Sharp2)  ',' num2str(SlewRate) ',' num2str(MedianWidth) ');' ] )
            disp(['peakfit(DataMatrix,'  num2str(X(round(xo))) ',' num2str(X(round(firstpoint))-X(round(lastpoint))) ')' ] )        
        case 111
            % When 'o' key is pressed, processed signal X,Y matrix is saved as in
            % mat file as the variable 'Output"
            Output=[X Y];
            uisave('Output');
        otherwise
            UnassignedKey=double(key)
            disp('Press k to print out list of keyboard commands')
    end % switch
end % if
% ----------------------------------------------------------------------
function [xo,dx]=panandzoom(X,xcenter,xrange)
 xo=val2ind(X,xcenter);
 hirange=val2ind(X,xcenter+xrange);
 lorange=val2ind(X,xcenter-xrange);
 dx=(hirange-lorange)./2;
 if xcenter<min(X),
      disp(['Lowest X value is ' num2str(min(X)) ]),
      xcenter=min(X)+xrange;
 end
 if xcenter>max(X),
       disp(['Highest X value is ' num2str(max(X)) ]),
       xcenter=max(X)-xrange;
  end
% ----------------------------------------------------------------------    
function [xx,yy]=RedrawSignal(x,y,xo,dx)
% Plots the entire signal (x,y) in the lower half of the plot window and an
% isolated segment (xx,yy) in the upper half, controlled by Pan and Zoom
% keys.
global X SmoothType SmoothWidth Sharpen Sharp1 Sharp2 ends
global DerivativeMode PeakLabels Overlay Report autozero MedianWidth
global SavedSignal GaussEstimate LorentzEstimate logymode SlewRate
Startx=round(xo-(dx/2));
Endx=abs(round(xo+(dx/2)-1));
if Endx>length(y),Endx=length(y);end
if Startx<1,Startx=1;end
PlotRange=Startx:Endx;
if (Endx-Startx)<2, PlotRange=xo:xo+2;end
xx=x(PlotRange);
yy=y(PlotRange); 
datasize=size(yy);if datasize(1)<datasize(2),yy=yy';end
datasize=size(xx);if datasize(1)<datasize(2),xx=xx';end
hold off
% auto-zero operation
lxx=length(xx);
bkgsize=5;
if autozero==1,
    X1=xx(1:round(lxx/bkgsize));
    X2=xx((lxx-round(lxx/bkgsize)):lxx);
    MX=[X1;X2];
    Y1=yy(1:round(length(xx)/bkgsize));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx); 
    MY=[Y1;Y2];
    bkgcoef=polyfit(MX,MY,1);  % Fit straight line to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if autozero
clf
% Plots isolated segment (xx,yy) in the upper half
figure(1);subplot(2,1,1);
if Overlay,
    hold on
    plot(xx,SavedSignal(PlotRange),'b:'); 
end % Overlay
plot(xx,yy,'b')
title('iSignal 2.5. Arrow keys to pan and zoom. Press K for keyboard commands')
lyy=min(yy);
uyy=max(yy)+(max(yy)-min(yy))/10;
if lyy<uyy;
   axis([x(Startx) x(Endx) lyy uyy ]);
end
center=x(round(xo));
hold on;plot([center center],[lyy uyy],'g-')
yyrange=max(yy)-min(yy);
if autozero==1,
  xlabel(sprintf('%s    y: %0.3g at %0.5g      P/P: %0.3g      Area: %0.3g    Std. Dev.: %0.2g ','Autozero ON ',y(round(xo)), center, yyrange, trapz(xx,yy), std(yy) ))
else
  xlabel(sprintf('%s    y: %0.3g at %0.5g      P/P: %0.2g      Area: %0.3g     Std. Dev.: %0.2g ','Autozero OFF',y(round(xo)), center, yyrange, trapz(xx,yy), std(yy) ))    
end
% Bottom half of the figure shows full signal
subplot(2,1,2);cla
if Overlay,
    hold on
    if logymode,
       semilogy(x,SavedSignal,'b:'); 
    else
       plot(x,SavedSignal,'b:'); 
    end % if logymode
end % Overlay
if logymode,
    semilogy(x,y,'b')  % Graph the signal with linear y axis on lower half
else
    plot(x,y,'b')  % Graph the signal with linear y axis on lower half
end % if logymode
axis([x(1) x(length(x)) min(min(y)) max(max(y))]); % Update plot
title('Smooth: S, A/Z   Derivatives: D   Sharpen: E,F/V G/B   Peak Meas: P   lin/log: H')
if Sharpen,
    xlabel([ num2str(SmoothWidth) ' point '  SmoothType ' smooth. Ends: ' num2str(ends) '   Der: ' num2str(DerivativeMode)  '   S1: ' num2str(round(10*Sharp1)/10)  '   S2: ' num2str(round(100*Sharp2)/100) '   Slew: ' num2str(SlewRate) '   Median: ' num2str(MedianWidth) ])
else
    xlabel([ num2str(SmoothWidth) ' point '  SmoothType ' smooth. Ends: ' num2str(ends) '    Der: ' num2str(DerivativeMode) '    Slew: ' num2str(SlewRate) '    Median: ' num2str(MedianWidth) ])    
end
if logymode,
    ylabel('Log y')
else
    ylabel('Linear y')
end
hold on
    % Mark the zoom range on the full signal with two magenta dotted vertical lines
checkzero=abs(y);
checkzero(~checkzero)=NaN; % Find smallest non-zero value
plot([min(xx) min(xx)],[min(checkzero) max(y)],'m--')
plot([max(xx) max(xx)],[min(checkzero) max(y)],'m--') 
plot([center center],[min(checkzero) max(y)],'g-')
if PeakLabels,
   % Compute limited range around center of zoom region
   FitLow=round(xo-(dx/10));
   if FitLow<1,FitLow=1;end
   FitHigh=abs(round(xo+(dx/10)-1));
   if FitHigh>length(X),FitHigh=length(X);end
   FitRange=FitLow:FitHigh;
   xxx=x(FitRange);
   yyy=y(FitRange);
   datasize=size(yyy);
   if datasize(1)<datasize(2),yyy=yyy';end
   datasize=size(xxx);
   if datasize(1)<datasize(2),xxx=xxx';end
   if autozero,
       yyy=yyy-polyval(bkgcoef,xxx);
   end
   [coef,S]=polyfit(xxx,yyy,2);  % Fit parabola to sub-group of points
   c1=coef(3);c2=coef(2);c3=coef(1);
   % Compute the correlation coefficient and R-Squared
   cc=corrcoef(polyval(coef,xxx),yyy);
   RSquared=cc(2).^2;
   % disp([c1 c2 c3 RSquared])  % for testing
   PeakX=-c2/(2*c3); % x-value of vertex
   PeakY=(c1-(c2*c2/(4*c3))); % y-value of vertex
   MeasuredWidth=norm(2.35703/(sqrt(2)*sqrt(-1*c3))); % Estimated Gaussian half-width
   % Label the peaks on the upper graph with number, position, height, and
   % width
   residual=yyy-polyval(coef,xxx);
   SNR=abs(PeakY./std(residual));
   if c2>0,
       % Fit parabola to log10 of sub-group
       [coef,S,MU]=polyfit(xxx,log(abs(yyy)),2);
       d1=coef(3);d2=coef(2);d3=coef(1);
       % Compute peak position and height of fitted parabola
       PeakX=-((MU(2).*d2/(2*d3))-MU(1));
       PeakY=exp(d1-d3*(d2/(2*d3))^2);
       MeasuredWidth=norm(MU(2).*2.35482/(sqrt(2)*sqrt(-1*d3)));
       % cc=corrcoef(polyval(coef,xxx),log(abs(yyy)));
       % RSquared=cc(2).^2;
       residual=yyy-PeakY*gaussian(xxx,PeakX,MeasuredWidth);
       SNR=abs(PeakY./std(residual));
   end
   subplot(2,1,1);
   hold on
   topaxis=axis;
   hpos=min(xx);
   yrange=topaxis(4)-topaxis(3);
   pos1=.1*yrange;
   pos2=.2*yrange;
   pos3=.3*yrange;
   pos4=.4*yrange;
   pos5=.5*yrange;
   pos6=.6*yrange;
   pos7=.7*yrange;
   text(hpos,topaxis(4)-pos1,[' Position=' num2str(PeakX)])
   text(hpos,topaxis(4)-pos2,[' Height=' num2str(PeakY)])
   text(hpos,topaxis(4)-pos3,[' Width=' num2str(MeasuredWidth)])
   text(hpos,topaxis(4)-pos4,[' Gaussian area=' num2str(1.0645*PeakY*MeasuredWidth) ])  
   area=trapz(xx,yy); % Compute the area
   text(hpos,topaxis(4)-pos5,[' Total area=' num2str(area) ])
   text(hpos,topaxis(4)-pos6,[' R2=' num2str(RSquared)])
   text(hpos,topaxis(4)-pos7,[' SNR=' num2str(round(10.*SNR)/10) ])
   if Report,
       % disp([PeakX PeakY MeasuredWidth 1.0645*PeakY*MeasuredWidth area SNR]);
       disp(sprintf('%0.5g       %0.5g       %0.5g       %0.5g       %0.5g       %0.3g',PeakX, PeakY, MeasuredWidth, 1.0645*PeakY*MeasuredWidth, area, SNR));
       Report=0;
   end % Report
   plotspace=linspace(min(xxx),max(xxx));
   if c2>0,
       plot(plotspace,PeakY.*exp(-((plotspace-PeakX)./(0.6005615.*MeasuredWidth)).^2),'r')
   else
       plot(plotspace,c3.*plotspace.^2+c2.*plotspace+c1,'r')
   end
   hold off
   xinterval=X(round(xo))-X(round(xo-1));
  if GaussEstimate,
    SmoothWidth=round(0.4*MeasuredWidth./xinterval);
    Sharp1=((MeasuredWidth/xinterval)^2)/25;
    Sharp2=((MeasuredWidth/xinterval)^4)/800;
    GaussEstimate=0;
  end % if GaussEstimate
  if LorentzEstimate,
    SmoothWidth=round(0.3*MeasuredWidth./xinterval);
    Sharp1=((MeasuredWidth/xinterval)^2)/8;
    Sharp2=((MeasuredWidth/xinterval)^4)/700;
    LorentzEstimate=0;
  end % if LorentzEstimate
end  % PeakLabels
 hold off
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
function Processed=ProcessSignal(x,y,DerivativeMode,w,type,ends,Sharpen,factor1,factor2,SlewRate,MedianWidth)
if SlewRate,
    for n=1:length(y)-1,
        if y(n+1)-y(n)>SlewRate,y(n+1)=y(n)+SlewRate;end
        if y(n+1)-y(n)<-SlewRate,y(n+1)=y(n)-SlewRate;end
    end
end % SlewRate
if MedianWidth,
    mY=y;
    for n=1:length(x)-(1+MedianWidth*2),
        mY(n+MedianWidth)=median(y((n):(n+1+MedianWidth*2)));
        y=mY;
    end
end  % MedianWidth
if type==0,w=1;end
if type==4,
    if w<2,w=3;end
    if DerivativeMode>4,
        if w<5,w=5;end
    end 
    % The polynomial order, 2+DerivativeMode, must be less than the
    % frame size, 2*w+1, and 2*w+1 must be odd.  
        Processed=savitzkyGolayFilt(y,2+DerivativeMode,DerivativeMode,2*w+1);
        if DerivativeMode==1,Processed=-Processed;end
        if DerivativeMode==3,Processed=-Processed;end;
else
    switch DerivativeMode
        case 0
            Processed=fastsmooth(y,w,type,ends);
        case 1
            Processed=fastsmooth(deriv(x,y),w,type,ends);
        case 2
            D2=fastsmooth(secderiv(x,y),w,1,ends);
            Processed=fastsmooth(D2,w,type,ends);
        case 3
            D2=fastsmooth(secderiv(x,y),w,1,ends);
            D3=fastsmooth(deriv(x,D2),w,1,ends);
            Processed=fastsmooth(D3,w,type,ends);
        case 4
            D2=fastsmooth(secderiv(x,y),w,1,ends);
            D4=fastsmooth(secderiv(x,D2),w,1,ends);
            D4=fastsmooth(D4,w,1,ends);
            
            Processed=fastsmooth(D4,w,type,ends);
        case 5
            D2=fastsmooth(secderiv(x,y),w,1,ends);
            D4=fastsmooth(secderiv(x,D2),w,1,ends);
            D4=fastsmooth(D4,w,1,ends);
            Processed=fastsmooth(deriv(x,D4),w,type,ends);
    end
end
if Sharpen,
    type=4; 
    if w<3;w=3;end
    Processed=enhance(x,Processed,factor1,factor2,w,type);
end
% following line for testing
% disp(['X,Y size = ' num2str(size(x)) ' , '  num2str(size(Processed)) ] )
Processed=reshape(Processed,size(x));
% datasize=size(Processed);if datasize(1)>datasize(2),Processed=Processed';end
% ----------------------------------------------------------------------
function Enhancedsignal=enhance(x,signal,factor1,factor2,SmoothWidth,type)
% Resolution enhancement function by even derivative method. The
% arguments factor1 and factor 2 are 2nd and 4th derivative weighting
% factors. Larger values of factor1 and factor2 will reduce the 
% peak width but will cause artifacts in the baseline near 
% the peak.  Adjust the factors for the the best compromise. 
% Use minimum smooth width needed to reduce excess noise. 
datasize=size(signal);
if datasize(1)>datasize(2),signal=signal';end
if type==4,
    Enhancedsignal = signal-factor1.*savitzkyGolayFilt(signal,4,2,2*SmoothWidth+1)+...
        factor2.*savitzkyGolayFilt(signal,6,4,2*SmoothWidth+1);
else
    d2=secderiv(x,signal);  % Computes second derivative
    d4=secderiv(x,d2);   % Computes fourth derivative
    Enhancedsignal = signal-factor1.*fastsmooth(d2,SmoothWidth,type)+...
        factor2.*fastsmooth(fastsmooth(fastsmooth(d4,SmoothWidth,2),SmoothWidth,2),SmoothWidth,2);
end
Enhancedsignal=Enhancedsignal';
% ----------------------------------------------------------------------
function d=secderiv(x,a)
% Second derivative of y with respect to x using 3-point central difference.
%  T. C. O'Haver, 2011.
n=length(a);
d=zeros(size(a));
for j = 2:n-2;
  x1=x(j-1);x2=x(j);x3=x(j+1);
  d(j)=((a(j+1)-a(j))./(x3-x2) - (a(j)-a(j-1))./(x2-x1))./((x3-x1)/2);
end
d(1)=d(2);
d(n)=d(n-1);
% ----------------------------------------------------------------------
function d=deriv(x,y)
% First derivative of y with respect to x using 2-point central difference.
%  T. C. O'Haver, 2011.
n=length(y);
d=zeros(size(y));
d(1)=(y(2)-y(1))./(x(2)-x(1));
d(n)=(y(n)-y(n-1))./(x(n)-x(n-1));
for j = 2:n-1;
  d(j)=(y(j+1)-y(j-1)) ./ (1.*(x(j+1)-x(j-1)));
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
    case 0
       SmoothY=sa(Y,w,ends);  
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
function sy=condense(y,n)
% Condense y by a factor of n, where n is a non-zero positive integer.
% Produces a shorter, approximate version of vector y, with each group 
% of n adjacent points in y replaced by its average. Use for reducing the 
% length and processing time of over-sampled signals or for preliminary 
% and exploratory analysis of very large signals to locate the interesting 
% bits, which can then be selected out of the full-length signal for 
% more precise analysis. For x,y data sets, use this function on both 
% independent variable x AND dependent variable y so that the features 
% of y will appear at the same x values.
% Example: condense([1 2 3 4 5 6 7 8 9 10 11 12],3) yields [2 5 8 11]
% condense([.9 1.1 .9 1 .9 1.1 .9 1 .9 1.1 .9 1],3) = [0.9667 1 0.9333 1]
% condense([0 0 0 0 0 0 0 1 1 1 1 1 1 1],2) = [0 0 0 .5 1 1 1]
n=round(n);
m=floor(length(y)/n);
if n > 1
    sy=mean(reshape(y(1:n*m),n,m));
else
    sy=y;
end
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)) .^2);
% ----------------------------------------------------------------------
function y=savitzkyGolayFilt(x,N,DN,F,W,DIM)
%savitzkyGolayFilt Savitzky-Golay Filtering.
%   savitzkyGolayFilt(X,N,DN,F) filters the signal X using a Savitzky-Golay 
%   (polynomial) filter.  The polynomial order, N, must be less than the
%   frame size, F, and F must be odd.  DN specifies the differentiation
%   order (DN=0 is smoothing). For a DN higher than zero, you'll have to
%   scale the output by 1/T^DN to acquire the DNth smoothed derivative of
%   input X, where T is the sampling interval. The length of the input X
%   must be >= F.  If X is a matrix, the filtering is done on the columns
%   of X.
%
%   Note that if the polynomial order N equals F-1, no smoothing
%   will occur.
%
%   savitzkyGolayFilt(X,N,DN,F,W) specifies a weighting vector W with
%   length F containing real, positive valued weights employed during the
%   least-squares minimization. If not specified, or if specified as
%   empty, W defaults to an identity matrix.
%
%   savitzkyGolayFilt(X,N,DN,F,[],DIM) or savitzkyGolayFilt(X,N,DN,F,W,DIM)
%   operates along the dimension DIM.
%   Copyright (c) 2011, Diederick
%   See also savitzkyGolay, FILTER, sgolayfilt

%   References:
%     [1] Sophocles J. Orfanidis, INTRODUCTION TO SIGNAL PROCESSING,
%              Prentice-Hall, 1995, Chapter 8.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.11.4.4 $  $Date: 2009/08/11 15:47:54 $

error(nargchk(4,6,nargin,'struct'));

% Check if the input arguments are valid
if round(F) ~= F, error(generatemsgid('MustBeInteger'),'Frame length must be an integer.'), end
if rem(F,2) ~= 1, error(generatemsgid('SignalErr'),'Frame length must be odd.'), end
if round(N) ~= N, error(generatemsgid('MustBeInteger'),'Polynomial order must be an integer.'), end
if N > F-1, error(generatemsgid('InvalidRange'),'The Polynomial order must be less than the frame length.'), end
if DN > N, error(generatemsgid('InvalidRange'),'The Differentiation order must be less than or equal to the Polynomial order.'), end

if nargin < 5 || isempty(W)
   % No weighting matrix, make W an identity
   W = ones(F,1);
else
   % Check for right length of W
   if length(W) ~= F, error(generatemsgid('InvalidDimensions'),'The weight vector must be of the same length as the frame length.'),end
   % Check to see if all elements are positive
   if min(W) <= 0, error(generatemsgid('InvalidRange'),'All the elements of the weight vector must be greater than zero.'), end
end

if nargin < 6, DIM = []; end

% Compute the projection matrix B
pp = fix(-F./2):fix(F./2);
B = savitzkyGolay(pp,N,DN,pp,W);

if ~isempty(DIM) && DIM > ndims(x)
	error(generatemsgid('InvalidDimensions'),'Dimension specified exceeds the dimensions of X.')
end

% Reshape X into the right dimension.
if isempty(DIM)
	% Work along the first non-singleton dimension
	[x, nshifts] = shiftdim(x);
else
	% Put DIM in the first dimension (this matches the order 
	% that the built-in filter function uses)
	perm = [DIM,1:DIM-1,DIM+1:ndims(x)];
	x = permute(x,perm);
end

if size(x,1) < F, error(generatemsgid('InvalidDimensions'),'The length of the input must be >= frame length.'), end

% Preallocate output
y = zeros(size(x));

% Compute the transient on (note, this is different than in sgolayfilt,
% they had an optimization leaving out some transposes that is only valid
% for DN==0)
y(1:(F+1)/2-1,:) = fliplr(B(:,(F-1)/2+2:end)).'*flipud(x(1:F,:));

% Compute the steady state output
ytemp = filter(B(:,(F-1)./2+1),1,x);
y((F+1)/2:end-(F+1)/2+1,:) = ytemp(F:end,:);

% Compute the transient off
y(end-(F+1)/2+2:end,:) = fliplr(B(:,1:(F-1)/2)).'*flipud(x(end-(F-1):end,:));

% Convert Y to the original shape of X
if isempty(DIM)
	y = shiftdim(y, -nshifts);
else
	y = ipermute(y,perm);
end
% ----------------------------------------------------------------------
function [fc, df] = savitzkyGolay(x,n,dn,x0,W,flag)
% Function:
%       Savitzky-Golay Smoothing and Differentiation Filter
%       Copyright (c) 2011, Diederick
%       The Savitzky-Golay smoothing/differentiation filter (i.e., the
%       polynomial smoothing/differentiation filter, or  the least-squares
%       smoothing/differentiation filters) optimally fit a set of data
%       points to polynomials of different degrees. 
%       See for details in Matlab Documents (help sgolay). The sgolay
%       function in Matlab can deal with only symmetrical and uniformly
%       spaced data of even number.
%       This function presented here is a general implement of the sgolay
%       function in Matlab. The Savitzky-Golay filter coefficients for even
%       number, nonsymmetrical and nonuniformly spaced data can be
%       obtained. And the filter coefficients for the initial point or the
%       end point can be obtained too. In addition, either numerical
%       results or symbolical results can be obtained. Lastly, this
%       function is faster than MATLAB's sgolay.
%
% Usage:
%       [fc,df] = savitzkyGolay(x,n,dn,x0,flag)
%   input:
%       x    = the original data point, e.g., -5:5 
%       n    = polynomial order
%       dn   = differentation order (0=smoothing),  default=0
%       x0   = estimation point, can be a vector    default=0
%       W    = weight vector, can be empty          
%              must have same length as x0          default=identity
%       flag = numerical(0) or symbolical(1),       default=0
%
%   output:
%       fc   = filter coefficients obtained (B output of sgolay).
%       df   = differentiation filters (G output of sgolay).
% Notes:
% 1.    x can be arbitrary, e.g., odd number or even number, symmetrical or
%       nonsymmetrical, uniformly spaced or nonuniformly spaced, etc.       
% 2.    x0 can be arbitrary, e.g., the initial point, the end point, etc.
% 3.    Either numerical results or symbolical results can be obtained.
% Example:
%       sgsdf([-3:3],2,0,0,[],0)
%       sgsdf([-3:3],2,0,0,[],1)
%       sgsdf([-3:3],2,0,-3,[],1)
%       sgsdf([-3:3],2,1,2,[],1)
%       sgsdf([-2:3],2,1,1/2,[],1)
%       sgsdf([-5:2:5],2,1,0,[],1)     
%       sgsdf([-1:1 2:2:8],2,0,0,[],1)
% Author:
%       Diederick C. Niehorster <dcniehorster@hku.hk> 2011-02-05
%       Department of Psychology, The University of Hong Kong
%
%       Originally based on
%       http://www.mathworks.in/matlabcentral/fileexchange/4038-savitzky-golay-smoothing-and-differentiation-filter
%       Allthough I have replaced almost all the code (partially based on
%       the comments on the FEX submission), increasing its compatibility
%       with MATLABs sgolay (now supports a weight matrix), its numerical
%       stability and it speed. Now, the help is pretty much all that
%       remains.
%       Jianwen Luo <luojw@bme.tsinghua.edu.cn, luojw@ieee.org> 2003-10-05
%       Department of Biomedical Engineering, Department of Electrical Engineering
%       Tsinghua University, Beijing 100084, P. R. China  
% Reference
%[1]A. Savitzky and M. J. E. Golay, "Smoothing and Differentiation of Data
%   by Simplified Least Squares Procedures," Analytical Chemistry, vol. 36,
%   pp. 1627-1639, 1964.
%[2]J. Steinier, Y. Termonia, and J. Deltour, "Comments on Smoothing and
%   Differentiation of Data by Simplified Least Square Procedures,"
%   Analytical Chemistry, vol. 44, pp. 1906-1909, 1972.
%[3]H. H. Madden, "Comments on Savitzky-Golay Convolution Method for
%   Least-Squares Fit Smoothing and Differentiation of Digital Data,"
%   Analytical Chemistry, vol. 50, pp. 1383-1386, 1978.
%[4]R. A. Leach, C. A. Carter, and J. M. Harris, "Least-Squares Polynomial
%   Filters for Initial Point and Slope Estimation," Analytical Chemistry,
%   vol. 56, pp. 2304-2307, 1984.
%[5]P. A. Baedecker, "Comments on Least-Square Polynomial Filters for
%   Initial Point and Slope Estimation," Analytical Chemistry, vol. 57, pp.
%   1477-1479, 1985.
%[6]P. A. Gorry, "General Least-Squares Smoothing and Differentiation by
%   the Convolution (Savitzky-Golay) Method," Analytical Chemistry, vol.
%   62, pp. 570-573, 1990.
%[7]Luo J W, Ying K, He P, Bai J. Properties of Savitzky-Golay Digital
%   Differentiators, Digital Signal Processing, 2005, 15(2): 122-136.
%
%See also:
%       sgolay, savitzkyGolayFilt

% Check if the input arguments are valid and apply defaults
error(nargchk(2,6,nargin,'struct'));

if round(n) ~= n, error(generatemsgid('MustBeInteger'),'Polynomial order (n) must be an integer.'), end
if round(dn) ~= dn, error(generatemsgid('MustBeInteger'),'Differentiation order (dn) must be an integer.'), end
if n > length(x)-1, error(generatemsgid('InvalidRange'),'The Polynomial Order must be less than the frame length.'), end
if dn > n, error(generatemsgid('InvalidRange'),'The Differentiation order must be less than or equal to the Polynomial order.'), end

% set defaults if needed
if nargin<6
    flag=false;
end
if nargin < 5 || isempty(W)
   % No weighting matrix, make W an identity
   W = eye(length(x0));
else
   % Check W is real.
   if ~isreal(W), error(generatemsgid('NotReal'),'The weight vector must be real.'),end
   % Check for right length of W
   if length(W) ~= length(x0), error(generatemsgid('InvalidDimensions'),'The weight vector must be of the same length as the frame length.'),end
   % Check to see if all elements are positive
   if min(W) <= 0, error(generatemsgid('InvalidRange'),'All the elements of the weight vector must be greater than zero.'), end
   % Diagonalize the vector to form the weighting matrix
   W = diag(W);
end
if nargin<4
    x0=0;
end
if nargin<3
    dn=0;
end

% prepare for symbolic output
if flag
    x=sym(x);
    x0=sym(x0);
end

Nx  = length(x);
x=x(:);
Nx0 = length(x0);
x0=x0(:);

if flag
    A=ones(length(x),1);
    for k=1:n
        A=[A x.^k];
    end
    df = inv(A'*A)*A';                          % backslash operator doesn't work as expected with symbolic inputs, but the "slowness and inaccuracy" of this method doesn't matter when doing the symbolic version
else
    df = cumprod([ones(Nx,1) x*ones(1,n)],2) \ eye(Nx);
end
df = df.';

hx = [(zeros(Nx0,dn)) ones(Nx0,1)*prod(1:dn)];  % order=0:dn-1,& dn,respectively
for k=1:n-dn                                    % order=dn+1:n=dn+k
    hx = [hx x0.^k*prod(dn+k:-1:k+1)];
end

% filter coeffs
fc = df*hx'*W;
