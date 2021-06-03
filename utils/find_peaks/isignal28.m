function pY=isignal(DataMatrix,xcenter,xrange,sm,sw,em,dm,rm,s1,s2,sr,mw)
% Y=isignal(DataMatrix,xcenter,xrange,SmoothMode,SmoothWidth,ends,...
% DerivativeMode,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth)
% Version 2.8. Modified frequency Spectrum function (Shift-S/Shift-A/Shift-X)
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
% sm,sw,em,dm,rm,s1,s2).  Last few versions:
% Version 2.2: Adds slew rate limit (~ key) and spike filter (M key).
% Version 2.3: Measures peak area two ways. Fixes several bugs related to input arguments.
% Version 2.4: Computes signal-to-noise ratio (SNR) in Peak mode.
% Version 2.5: Bug fixes
% Version 2.6: Added peakfit function (Shift-F)
% Version 2.7: Added frequenecy Spectrum function (Shift-S/Shift-A/Shift-Z)
% By T. C. O'Haver (toh@umd.edu); Savitzky-Golay smooth code by Diederick.  
% See http://terpconnect.umd.edu/~toh/spectrum/iSignal.html
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
% Example 10: Weak peak at x=128 on a smooth, curved background. 
% Try second derivative + smoothing
% x=1:.1:256;
% y=gaussian(x,-100,300)+.02.*gaussian(x,128,30)+0.001.*randn(size(x));
% isignal(x,y);
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
%  Cycle derivative orders.....D/Shift-D Increase/Decrease derivative order
%  Toggle peak sharpening......E (0=OFF 1=ON)
%  Sharpening for Gaussian.....Y  Set sharpen settings for Gaussian
%  Sharpening for Lorentzian...U  Set sharpen settings for Lorentzian
%  Adjust sharp1...............F,V  F=>sharper, V=>less sharpening
%  Adjust sharp2...............G,B  G=>sharper, B=>less sharpening
%  Slew rate limit (0=OFF).....`  Largest allowed change between points
%  Spike filter width (O=OFF)..m  Spike filter eliminates sharp spikes
%  Toggle peak parabola........P  fits parabola to center, labels vertex
%  Fits peak in upper window...Shift-F (Asks for shape, number of peaks, etc)
%  Spectrum mode on/off........Shift-S (Shift-A and Shift-X to change axes)
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
global X Y xo dx DerivativeMode Sharpen Sharp1 Sharp2 SmoothWidth SlewRate MedianWidth
global SmoothType ends SmoothMode SavedSignal SavedXvalues PeakLabels Report autozero
global plotmode xmode SpectrumMode
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
plotmode=2; % Frequency spectrum initially in semilog y mode.
NumPeaksUW=1;
xmode=0;% Frequency spectrum initially in frequency mode.
SpectrumMode=0; % Frequency spectrum initially off. 
            
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
global X Y xx yy xo dx SmoothMode SmoothWidth DerivativeMode SlewRate extra SpectrumMode
global Sharpen Sharp1 Sharp2 SavedSignal SavedXvalues SavedBackground plotmode
global SmoothType ends PeakLabels Overlay Report GaussEstimate MedianWidth xmode
global LorentzEstimate logymode autozero Shape NumTrials NumPeaksUW fixedparameters
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
        case 68
            % When 'Shift-d' key is pressed, cycles through DerivativeModes 0,1,2,3,4,5->0
            % if length(Y)>10000,disp('Warning: Derivatives can be slow for for signal lengths above 10,000 points'),end
            DerivativeMode=DerivativeMode-1;
            if DerivativeMode==-1,DerivativeMode=5; end
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
%       case 95  % FUTURE ADDITION?
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
            disp('KEYBOARD CONTROLS, version 2.7:')
            disp(' Pan signal left and right...Coarse pan: < and >')
            disp('                             Fine pan: left and right cursor arrows')
            disp('                             Nudge: [ and ] ')
            disp(' Zoom in and out.............Coarse zoom: / and "  ')
            disp('                             Fine zoom: up and down cursor arrows')
            disp(' Resets pan and zoom.........ESC')
            disp(' Adjust smooth width.........A,Z (A=>more, Z=>less) ')
            disp(' Cycle smooth types..........S (No, Rectanglular, Triangle, Gaussian, Savitzky-Golay)')
            disp(' Toggle smooth ends..........X (0=ends zeroed  1=ends smoothed (slower)')
            disp(' Cycle derivative orders.....D/Shift-D Increase/Decrease derivative order')
            disp(' Toggle peak sharpening......E (0=OFF 1=ON)')
            disp(' Sharpening for Gaussian.....Y  Set sharpen settings for Gaussian')
            disp(' Sharpening for Lorentzian...U  Set sharpen settings for Lorentzian')
            disp(' Adjust sharp1...............F,V  F=>sharper, V=>less sharpening')
            disp(' Adjust sharp2   ............G,B  G=>sharper, B=>less sharpening')
            disp(' Slew rate limit (0=OFF).....~  Largest allowed y change between points')
            disp(' Spike filter width (O=OFF)..m  spike filter eliminates sharp spikes')
            disp(' Toggle peak parabola........P  fits parabola to center, labels vertex')
            disp(' Fits peak in upper window...Shift-F (Asks for shape, number of peaks, etc)')
            disp(' Spectrum mode on/off........Shift-S (Shift-A and Shift-X to change axes)')            
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
        case 70
            % When 'Shift-F' is pressed, applies peakfit function only to
            %  peaks in the upper window.
            Startx=round(xo-(dx/2));
            Endx=abs(round(xo+(dx/2)-1));
            if Endx>length(Y),Endx=length(Y);end
            if Startx<1,Startx=1;end
            PlotRange=Startx:Endx;
            if (Endx-Startx)<2, PlotRange=xo:xo+2;end
            xx=X(PlotRange);
            yy=Y(PlotRange);
            disp(' ')
            disp('Select the peak shape of the model (type 1-17 and press Enter key):')
            disp('Gaussian: y=exp(-((x-pos)./(0.6005615.*width)) .^2)')
            disp('  Gaussians with independent positions and widths : 1 (default)')
            disp('  Exponentionally-broadened Gaussian : 5 ')
            disp('  Exponentionally-broadened equal-width Gaussian : 8')
            disp('  Equal-width Gaussians : 6')
            disp('  Fixed-width Gaussians : 11')
            disp('  Fixed-position Gaussians : 16 ')
            disp('  Asymmetrical Gaussians (unequal half-widths on both sides) : 14.')
            disp('Lorentzian: y=ones(size(x))./(1+((x-pos)./(0.5.*width)).^2)')
            disp('  Lorentzians with independent positions and widths : 2')
            disp('  Equal-width Lorentzians : 7')
            disp('  Fixed-width Lorentzians : 12')
            disp('  Fixed-position Lorentzians : 17')
            disp('  Asymmetrical Lorentzians (unequal half-widths on both sides) : 15')
            disp('Blended sum of Gaussian and Lorentzian functions : 13')
            disp('Logistic: n=exp(-((x-pos)/(.477.*wid)).^2); y=(2.*n)./(1+n) : 3  ')
            disp('Pearson: y=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m : 4')
            disp('Exponential pulse: y=exp(-tau1.*x).*(1-exp(-tau2.*x)) : 9')
            disp('Sigmoid: y=1/2 + 1/2* erf((x-tau1)/sqrt(2*tau2)) : 10')
            disp(' ')
            Shapeinput=input('Peak shape number: ');
            if isempty(Shape),
            else
                Shape=Shapeinput;
            end
            if Shape>17, Shape=17;end
            if Shape<1, Shape=1;end
            if isempty(Shape),Shape=1;end
            switch Shape
                case 1
                    ShapeString='Gaussian';
                case 2
                    ShapeString='Lorentzian';
                case 3
                    ShapeString='logistic';
                    inputextra=input('Shape factor: ');
                    if isempty(inputextra),
                    else
                        extra=inputextra;
                    end
                case 4
                    ShapeString='Pearson';
                case 5
                    ShapeString='ExpGaussian';
                    inputextra=input('Exponentional factor: ');
                    if isempty(inputextra),
                    else
                        extra=inputextra;
                    end
                case 6
                    ShapeString='Equal-width Gaussian';
                case 7
                    ShapeString='Equal-width Lorentzian';
                case 8
                    ShapeString='Equal-width ExpGauss.';
                    inputextra=input('Exponentional factor: ');
                    if isempty(inputextra),
                    else
                        extra=inputextra;
                    end
                case 9
                    ShapeString='Exponental pulse';
                case 10
                    ShapeString='Sigmoid';
                case 11
                    ShapeString='Fixed-width Gaussian';
                    inputwidth=input('Peak width: ');
                    if isempty(inputwidth),
                    else
                        fixedparameters=inputwidth;
                    end
                case 12
                    ShapeString='Fixed-width Lorentzian';
                    inputwidth=input('Peak width: ');
                    if isempty(inputwidth),
                    else
                        fixedparameters=inputwidth;
                    end
                case 13
                    ShapeString='Gauss/Lorentz blend';
                    inputextra=input('Percent Gaussian: ');
                    if isempty(inputextra),
                    else
                        extra=inputextra;
                    end
                case 14
                    ShapeString='bifurcated Gaussian';
                    inputextra=input('Asymmetry factor: ');
                    if isempty(inputextra),
                    else
                        extra=inputextra;
                    end
                case 15
                    ShapeString='bifurcated Lorentzian';
                    inputextra=input('Asymmetry factor: ');
                    if isempty(inputextra),
                    else
                        extra=inputextra;
                    end
                case 16
                    ShapeString='Fixed-position Gaussians';
                    inputpositions=input('Peak positions as a vector, e.g. [200 400 600]: ');
                    if isempty(inputpositions),
                    else
                        fixedparameters=inputpositions;
                    end
                case 17
                    ShapeString='Fixed-position Lorentzians';
                    inputpositions=input('Peak positions as a vector, e.g. [200 400 600]: ');
                    if isempty(inputpositions),
                    else
                        fixedparameters=inputpositions;
                    end
                otherwise
                    ShapeString='';
            end % switch Shape''
            if Shape==11||Shape==12||Shape==16||Shape==17, % Only for fixed shapes
                NumPeaksUW=length(fixedparameters);
            else
                inputNumPeaks=input('Number of peaks: ');
                if isempty(inputNumPeaks),
                else
                    NumPeaksUW=inputNumPeaks;
                end
            end
            inputNumTrials=input('Number of Trial fits: ');
            if isempty(inputNumTrials),
            else
                if isnumeric(inputNumTrials),
                    NumTrials=inputNumTrials;
                else
                    NumTrials=1;
                end
            end
            X1=min(xx);
            X2=max(xx);
            lyy=min(yy);
            uyy=max(yy)+(max(yy)-min(yy))/10;
            n=X2-X1;
            width=n/(5*NumPeaksUW);
            startvector=[];
            startpos=[n/(NumPeaksUW+1):n/(NumPeaksUW+1):n-(n/(NumPeaksUW+1))]+X1;
            for marker=1:NumPeaksUW,
                markx=startpos(marker);
                startvector=[startvector markx width];
                plot([markx markx],[lyy uyy],'m--')
            end % for marker
            disp(['Least-squares fit of selected peaks to ' ShapeString ' peak model using the peakfit function:' ])
            figure(2)
            [FitResults,MeanFitError]=peakfit([xx,yy],0,0,NumPeaksUW,Shape,extra,NumTrials,startvector,autozero,fixedparameters);
            disp(['Fitting Error ' num2str(MeanFitError) '%'])
            disp('          Peak#     Position     Height      Width         Area  ')
            % for peak=1:NumPeaksUW,FitResults(peak,1)=PUW(peak,1);end
            disp(FitResults(:,1:5))
            disp('Peakfit plot shown in Figure 2')
            figure(1)
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 83 % If Shift-S is pressed, plots frequency spectrum in the lower window
            if SpectrumMode==1, 
                SpectrumMode=0;
                [xx,yy]=RedrawSignal(X,Y,xo,dx);
            else
                SpectrumMode=1;
                % Plot the power spectrum  in the lower half of
                % the window.
                PlotFrequencySpectrum(X,Y,xo,dx,plotmode,xmode)
                subplot(2,1,1)
                title('iSignal 2.8 Frequency Spectrum Mode (Press Shift-S again to cancel')
            end  
        case 65 % If Shift-A is pressed, changes plot mode for Spectrum
            plotmode=plotmode+1;
            if plotmode==5;plotmode=1;end
            PlotFrequencySpectrum(X,Y,xo,dx,plotmode,xmode)
            subplot(2,1,1)
            title('iSignal 2.8 Frequency Spectrum Mode (Press Shift-S again to cancel')
         case 88 % If Shift-X is pressed,changes xmode for Spectrum
            xmode=xmode+1;
            if xmode==2;xmode=0;end
            PlotFrequencySpectrum(X,Y,xo,dx,plotmode,xmode)
            subplot(2,1,1)
           title('iSignal 2.8 Frequency Spectrum Mode (Press Shift-S again to cancel')
        otherwise
            UnassignedKey=double(key)
            disp('Press k to print out list of keyboard commands')
    end % switch double(key),
end % if  isscalar(key),
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
global X SmoothType SmoothWidth Sharpen Sharp1 Sharp2 ends SpectrumMode
global DerivativeMode PeakLabels Overlay Report autozero MedianWidth xmode
global SavedSignal GaussEstimate LorentzEstimate logymode SlewRate plotmode 
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
title('iSignal 2.7. Arrow keys to pan and zoom. Press K for keyboard commands')
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
if SpectrumMode,
    PlotFrequencySpectrum(x,y,xo,dx,plotmode,xmode)
else
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
end % if SpectrumMode
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
% ----------------------------------------------------------------------
function PlotFrequencySpectrum(X,Y,xo,dx,plotmode,XMODE)
global SmoothType SmoothWidth Sharpen Sharp1 Sharp2 ends 
global DerivativeMode  MedianWidth SlewRate 
Startx=round(xo-(dx/2));
Endx=abs(round(xo+(dx/2)-1));
if Endx>length(Y),Endx=length(Y);end
if Startx<1,Startx=1;end
PlotRange=Startx:Endx;
if (Endx-Startx)<2, PlotRange=xo:xo+2;end
xx=X(PlotRange);
yy=Y(PlotRange);
subplot(2,1,1)
title('iSignal 2.8 Frequency Spectrum Mode (Press Shift-S again to cancel)')
xlabel('Press Shift-A to cycle through spectrum log/linear plot modes ')
subplot(2,1,2)
fy=fft(yy);
sy=fy .* conj(fy); % Compute power spectrum
plotrange=2:length(fy)/2;
if XMODE,
    f=range(xx)./(plotrange-1);
else
    f=((plotrange-1)./range(xx));
end
realsy=real(sy(plotrange));
maxpower=max(realsy);
maxx=val2ind(realsy,maxpower);
maxf=f(maxx);
hold off
switch plotmode,
    case 1
        plot(f,realsy,'r.-')
        ylabel('Linear y')
    case 2
        semilogx(f,realsy,'r.-')
        ylabel('Linear y')
    case 3
        semilogy(f,realsy,'r.-')
        ylabel('Log y')
    case 4
        loglog(f,realsy,'r.-')
        ylabel('Log y')
    otherwise,
end
spectrumaxis=axis;
hpos=min(realsy);
text(spectrumaxis(1),0.5.*spectrumaxis(4),['Peak ' num2str(maxf) ' at harmonic #' num2str(maxx) ])
 if Sharpen,
        xlabel([ num2str(SmoothWidth) ' point '  SmoothType ' smooth. Ends: ' num2str(ends) '   Der: ' num2str(DerivativeMode)  '   S1: ' num2str(round(10*Sharp1)/10)  '   S2: ' num2str(round(100*Sharp2)/100) '   Slew: ' num2str(SlewRate) '   Median: ' num2str(MedianWidth) ])
    else
        xlabel([ num2str(SmoothWidth) ' point '  SmoothType ' smooth. Ends: ' num2str(ends) '    Der: ' num2str(DerivativeMode) '    Slew: ' num2str(SlewRate) '    Median: ' num2str(MedianWidth) ])
 end
if XMODE,
    title('x=Time. Press Shift-X to change to frequency.')
else
    title('x=Frequency (e.g. 1/time). Press Shift-X to change to time.')
end
% ----------------------------------------------------------------------
function [FitResults,LowestError,BestStart,xi,yi,BootResults]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,fixedparameters,plots)
global AA xxx PEAKHEIGHTS fixedparameters
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
xoffset=center-window/2;
xoffset=0;
n1=val2ind(X,center-window/2);
n2=val2ind(X,center+window/2);
if window==0,n1=1;n2=length(X);end
xx=X(n1:n2)-xoffset;
yy=Y(n1:n2);
ShapeString='Gaussian';

% Define values of any missing arguments
switch nargin
    case 1
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        xx=X;yy=Y;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=1;
        plots=1;
    case 2
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        xx=signal;yy=center;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=1;
        plots=1;
    case 3
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=1;
        fixedparameters=0;
        plots=1;
    case 4
        peakshape=1;
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=1;
        fixedparameters=0;
        plots=1;
    case 5
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=1;
        fixedparameters=0;
        plots=1;
    case 6
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=1;
        fixedparameters=0;
        plots=1;
    case 7
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=1;
        fixedparameters=0;
        plots=1;
    case 8
        AUTOZERO=1;
        fixedparameters=0;
        plots=1;
    case 9
        fixedparameters=0;
        plots=1;
    case 10
        plots=1;
    case 11
        
    otherwise
end % switch nargin

% Default values for placeholder zeros
if NumTrials==0;NumTrials=1;end
if peakshape==0;peakshape=1;end
if NumPeaks==0;NumPeaks=1;end
if start==0;start=calcstart(xx,NumPeaks,xoffset);end
% if FIXEDWIDTH==0, FIXEDWIDTH=length(xx)/10;end
% if peakshape==16;FIXEDPOSITIONS=fixedparameters;end

% Remove linear baseline from data segment if AUTOZERO==1
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
% for peaks=1:NumPeaks,
%      peakindex=2*peaks-1;
%      newstart(peakindex)=start(peakindex)-xoffset;
% end

% Assign ShapStrings
switch peakshape
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
        ShapeString='Sigmoid';
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
    otherwise
end % switch peakshape
  
% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.001,'Display','off' );
LowestError=1000; % or any big number greater than largest error expected
FitParameters=zeros(1,NumPeaks.*2); 
BestStart=zeros(1,NumPeaks.*2); 
height=zeros(1,NumPeaks); 
bestmodel=zeros(size(yy));
for k=1:NumTrials, 
    % disp(['Trial number ' num2str(k) ] ) % optionally prints the current trial number as progress indicator
  switch peakshape
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
          TrialParameters=fminsearch(@(lambda)(fitsigmoid(lambda,xx,yy)),newstart,options);
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
              fixedstart(pc)=fixedstart(pc)+.1*randn().*fixedstart(pc);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart,options);
      case 17
           fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
              fixedstart(pc)=fixedstart(pc)+.1*randn().*fixedstart(pc);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart,options);
      otherwise
  end % switch peakshape
  
% Construct model from Trial parameters
A=zeros(NumPeaks,n);
for m=1:NumPeaks,
    switch peakshape
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
            A(m,:)=sigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 11
            A(m,:)=gaussian(xx,TrialParameters(m),fixedparameters);
        case 12
            A(m,:)=lorentzian(xx,TrialParameters(m),fixedparameters);
        case 13
            A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 14
            A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 15
            A(m,:)=BiLorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
        case 16
            A(m,:)=gaussian(xx,fixedparameters(m),TrialParameters(m));
        case 17
            A(m,:)=lorentzian(xx,fixedparameters(m),TrialParameters(m));
        otherwise
    end % switch
end % for

% Multiplies each row by the corresponding amplitude and adds them up
model=PEAKHEIGHTS'*A;

% Compare trial model to data segment and compute the fit error
  MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
  % Take only the single fit that has the lowest MeanFitError
  if MeanFitError<LowestError, 
      if min(PEAKHEIGHTS)>0,  % Consider only fits with positive peak heights
        LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
        FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
        BestStart=newstart; % Assign BestStart to the start with the lowest MeanFitError
        height=PEAKHEIGHTS; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
        bestmodel=model; % Assign bestmodel to the model with the lowest MeanFitError
      end % if min(PEAKHEIGHTS)>0
  end % if MeanFitError<LowestError
end % for k (NumTrials)
%
% Construct model from best-fit parameters
AA=zeros(NumPeaks,200);
xxx=linspace(min(xx),max(xx),200);
for m=1:NumPeaks,
   switch peakshape
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
        AA(m,:)=sigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m));   
    case 11
        AA(m,:)=gaussian(xxx,FitParameters(m),fixedparameters);
    case 12
        AA(m,:)=lorentzian(xxx,FitParameters(m),fixedparameters);
    case 13
        AA(m,:)=GL(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 14
        AA(m,:)=BiGaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 15
        AA(m,:)=BiLorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 16
        AA(m,:)=gaussian(xxx,fixedparameters(m),FitParameters(m));
    case 17
        AA(m,:)=lorentzian(xxx,fixedparameters(m),FitParameters(m));
       otherwise
  end % switch
end % for

% Multiplies each row by the corresponding amplitude and adds them up
heightsize=size(height');
AAsize=size(AA);
if heightsize(2)==AAsize(1),
   mmodel=height'*AA;
else
    mmodel=height*AA;
end
% Top half of the figure shows original signal and the fitted model.
if plots,
    subplot(2,1,1);plot(xx+xoffset,yy,'b.'); % Plot the original signal in blue dots
    hold on
end
for m=1:NumPeaks,
    if plots, plot(xxx+xoffset,height(m)*AA(m,:),'g'),end  % Plot the individual component peaks in green lines
    area(m)=trapz(xxx+xoffset,height(m)*AA(m,:)); % Compute the area of each component peak using trapezoidal method
    yi(m,:)=height(m)*AA(m,:); % (NEW) Place y values of individual model peaks into matrix yi
end
xi=xxx+xoffset; % (NEW) Place the x-values of the individual model peaks into xi

if plots,
    % Mark starting peak positions with vertical dashed lines
    if peakshape==16||peakshape==17
    else
        for marker=1:NumPeaks,
            markx=BestStart((2*marker)-1);
            subplot(2,1,1);plot([markx+xoffset markx+xoffset],[0 max(yy)],'m--')
        end % for
    end % if peakshape
    plot(xxx+xoffset,mmodel,'r');  % Plot the total model (sum of component peaks) in red lines
    hold off;
    axis([min(xx)+xoffset max(xx)+xoffset min(yy) max(yy)]);
    switch AUTOZERO,
        case 0
            title('Peakfit 3.6 Autozero OFF.')
        case 1
            title('Peakfit 3.6 Linear autozero.')
        case 2
            title('Peakfit 3.6 Quadratic autozero.')
    end
    if peakshape==4||peakshape==5||peakshape==8||peakshape==13||peakshape==14||peakshape==15, % Shapes with Extra factor
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '%    Extra = ' num2str(extra) ] )
    else
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '%' ] )
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
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models only
            FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12, % Fixed-width shapes only
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) fixedparameters area(m)]];
            else
                if peakshape==16||peakshape==17, % Fixed-position shapes only
                    FitResults=[round(m) fixedparameters(m) height(m) FitParameters(m) area(m)];
                else
                    FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)];
                end
            end
        end % if peakshape
    else
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models only
            FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12, % Fixed-width shapes only
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) fixedparameters area(m)]];
            else
                if peakshape==16||peakshape==17, % Fixed-position shapes only
                    FitResults=[FitResults ; [round(m) fixedparameters(m) height(m) FitParameters(m) area(m)]];
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
    if peakshape==9||peakshape==10,  % Pulse and sigmoid shapes only
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
end % if plots

if NumArgOut==6,
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
global PEAKHEIGHTS fixedparameters
format short g
format compact
warning off all
% FIXEDWIDTH=fixedparameters;
xoffset=0;
if start==0;start=calcstart(xx,NumPeaks,xoffset);end
PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
newstart=start;

for peaks=1:NumPeaks,
     peakindex=2*peaks-1;
     newstart(peakindex)=start(peakindex)-xoffset;
end

% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.00001,'Display','off' );
LowestError=1000; % or any big number greater than largest error expected
FitParameters=zeros(1,NumPeaks.*2); 
BestStart=zeros(1,NumPeaks.*2); 
height=zeros(1,NumPeaks); 
bestmodel=zeros(size(yy));
for k=1:NumTrials, 
  switch peakshape
    case 1
        TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart);
    case 2
        TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart);
    case 3
        TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart);
    case 4
        TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart);
    case 5
        zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
        zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
        TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart);
    case 6
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart);
    case 7
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitlorentziancw(lambda,xx,yy)),cwnewstart);
    case 8
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart);
      case 9
          TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart);
      case 10
          TrialParameters=fminsearch(@(lambda)(fitsigmoid(lambda,xx,yy)),newstart);
      case 11
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFWGaussian(lambda,xx,yy)),fixedstart);
      case 12
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFWLorentzian(lambda,xx,yy)),fixedstart);
      case 13
          TrialParameters=fminsearch(@(lambda)(fitGL(lambda,xx,yy,extra)),newstart);
      case 14
          TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart);
      case 15
          TrialParameters=fminsearch(@(lambda)(fitBiLorentzian(lambda,xx,yy,extra)),newstart);
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
      otherwise
  end % switch peakshape
  
% Construct model from Trial parameters
A=zeros(NumPeaks,n);
for m=1:NumPeaks,
   switch peakshape
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
        A(m,:)=sigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m)); 
    case 11
        A(m,:)=gaussian(xx,TrialParameters(m),fixedparameters);
    case 12
        A(m,:)=lorentzian(xx,TrialParameters(m),fixedparameters); 
    case 13
        A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
    case 14
        A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);       
    case 15
        A(m,:)=BiLorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);       
    case 16
        A(m,:)=gaussian(xx,fixedparameters(m),TrialParameters(m));
    case 17
        A(m,:)=lorentzian(xx,fixedparameters(m),TrialParameters(m));

   end % switch
    for parameter=1:2:2*NumPeaks,
        newstart(parameter)=newstart(parameter)*(1+randn/50);
        newstart(parameter+1)=newstart(parameter+1)*(1+randn/10);
    end
end % for
% Multiplies each row by the corresponding amplitude and adds them up
model=PEAKHEIGHTS'*A;

% Compare trial model to data segment and compute the fit error
  MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
  % Take only the single fit that has the lowest MeanFitError
  if MeanFitError<LowestError, 
      if min(PEAKHEIGHTS)>0,  % Consider only fits with positive peak heights
        LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
        FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
%         BestStart=newstart; % Assign BestStart to the start with the lowest MeanFitError
        height=PEAKHEIGHTS; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
%         bestmodel=model; % Assign bestmodel to the model with the lowest MeanFitError
      end % if min(PEAKHEIGHTS)>0
  end % if MeanFitError<LowestError
end % for k (NumTrials)

for m=1:NumPeaks,
    area(m)=trapz(xx+xoffset,height(m)*A(m,:)); % Compute the area of each component peak using trapezoidal method
end

for m=1:NumPeaks,
    if m==1,
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models
            FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12,  % Fixed-width shapes only
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) fixedparameters area(m)]];
            else
                FitResults=[[round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]];
            end
        end % if peakshape
    else
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models
            FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12, % Fixed-width shapes only
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) fixedparameters area(m)]];
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
function err = fitgaussian(lambda,t,y)
% Fitting function for a Gaussian band signal.
global PEAKHEIGHTS
numpeaks=round(length(lambda)/2);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(2*j-1),lambda(2*j))';
end 
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitewgaussian(lambda,t,y)
% Fitting function for a Gaussian band signal with equal peak widths.
global PEAKHEIGHTS
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),lambda(numpeaks+1))';
end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFWGaussian(lambda,t,y)
%	Fitting function for a fixed width Gaussian
global PEAKHEIGHTS fixedparameters
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),fixedparameters)';
end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFPGaussian(lambda,t,y)
%	Fitting function for fixed-position Gaussians
global PEAKHEIGHTS fixedparameters
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,fixedparameters(j), lambda(j))';
end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFPLorentzian(lambda,t,y)
%	Fitting function for fixed-position Lorentzians
global PEAKHEIGHTS fixedparameters
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,fixedparameters(j), lambda(j))';
end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFWLorentzian(lambda,t,y)
%	Fitting function for fixed width Gaussians
global PEAKHEIGHTS fixedparameters
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),fixedparameters)';
end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitlorentziancw(lambda,t,y)
% Fitting function for a Lorentzian band signal with equal peak widths.
global PEAKHEIGHTS
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),lambda(numpeaks+1))';
end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitlorentzian(lambda,t,y)
%	Fitting function for single lorentzian, lambda(1)=position, lambda(2)=width
%	Fitgauss assumes a lorentzian function 
global PEAKHEIGHTS
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = lorentzian(t,lambda(2*j-1),lambda(2*j))';
end
PEAKHEIGHTS = A\y';
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
%	Fitting function for logistic, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fitlogistic assumes a logistic function 
%  T. C. O'Haver, May 2006
global PEAKHEIGHTS
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = logistic(t,lambda(2*j-1),lambda(2*j))';
end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = logistic(x,pos,wid)
% logistic function.  pos=position; wid=half-width (both scalar)
% logistic(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991 
n = exp(-((x-pos)/(.477.*wid)) .^2);
g = (2.*n)./(1+n);
% ----------------------------------------------------------------------
function g = lognormal(x,pos,wid)
% lognormal function.  pos=position; wid=half-width (both scalar)
% lognormal(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991  
g = exp(-(log(x/pos)/(0.01.*wid)) .^2);
% ----------------------------------------------------------------------
function err = fitpearson(lambda,t,y,shapeconstant)
%   Fitting functions for a Pearson 7 band signal.
% T. C. O'Haver (toh@umd.edu),   Version 1.3, October 23, 2006.
global PEAKHEIGHTS
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = pearson(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
PEAKHEIGHTS = A\y';
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
global PEAKHEIGHTS
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = expgaussian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitexpewgaussian(lambda,t,y,timeconstant)
% Fitting function for exponentially-broadened Gaussian bands with equal peak widths.
global PEAKHEIGHTS
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = expgaussian(t,lambda(j),lambda(numpeaks+1),timeconstant);
end
PEAKHEIGHTS = abs(A\y');
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
function yb = ExpBroaden(y,t)
% ExpBroaden(y,t) zero pads y and convolutes result by an exponential decay
% of time constant t by multiplying Fourier transforms and inverse
% transforming the result.
hly=round(length(y)./2);
ey=[zeros(1,hly)';y;zeros(1,hly)'];
fy=fft(ey);
a=exp(-(1:length(fy))./t);
fa=fft(a);
fy1=fy.*fa';
ybz=real(ifft(fy1))./sum(a);
yb=ybz(hly+2:length(ybz)-hly+1);
% ----------------------------------------------------------------------
function err = fitexppulse(tau,x,y)
% Iterative fit of the sum of exponental pulses
% of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
global PEAKHEIGHTS
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = exppulse(x,tau(2*j-1),tau(2*j));
end
PEAKHEIGHTS =abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = exppulse(x,t1,t2)
% Exponential pulse of the form 
% Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
e=(x-t1)./t2;
p = 4*exp(-e).*(1-exp(-e));
p=p .* [p>0];
g = p';
% ----------------------------------------------------------------------
function err = fitsigmoid(tau,x,y)
% Fitting function for iterative fit to the sum of
% sigmiods of the form Height./(1 + exp((t1 - t)/t2))
global PEAKHEIGHTS
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = sigmoid(x,tau(2*j-1),tau(2*j));
end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g=sigmoid(x,t1,t2)
% g=1./(1 + exp((t1 - x)./t2))'; % old version of sigmoid
g=1/2 + 1/2* erf(real((x-t1)/sqrt(2*t2))); % Modified sigmoid
% Bifurcated sigmoid
% lx=length(x);
% hx=val2ind(x,t1);
% g(1:hx)=1./(1 + exp((t1 - x(1:hx))./t2));
% g(hx+1:lx)=1./(1 + exp((t1 - x(hx+1:lx))./(1.3*t2)));
% ----------------------------------------------------------------------
function err = fitGL(lambda,t,y,shapeconstant)
%   Fitting functions for Gaussian/Lorentzian blend.
% T. C. O'Haver (toh@umd.edu), 2012.
global PEAKHEIGHTS
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = GL(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = GL(x,pos,wid,m)
% Gaussian/Lorentzian blend. m = percent Gaussian character
% pos=position; wid=half-width
% m = percent Gaussian character.
%  T. C. O'Haver, 2012
g=((m/100)*gaussian(x,pos,wid)+(1-(m/100))*lorentzian(x,pos,wid))/2;
% ----------------------------------------------------------------------
function err = fitBiGaussian(lambda,t,y,shapeconstant)
%   Fitting functions for BiGaussian.
% T. C. O'Haver (toh@umd.edu),  2012.
global PEAKHEIGHTS
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = BiGaussian(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
PEAKHEIGHTS = A\y';
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
global PEAKHEIGHTS
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = BiLorentzian(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
PEAKHEIGHTS = A\y';
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