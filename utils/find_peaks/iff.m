function iff(x,y)
% Interactive Fourier Filter function for data in arguments x,y,
% with sliders for interactive control of the center frequency, 
% width, and shape of the filter. There are four operating modes:
% mode=0 for band-pass filter, linear plot; 
% mode=1 for  band-pass filter, semilogx plot; 
% mode=2 for band-reject (notch) filter, linear plot;
% mode=3 for band-reject (notch) filter, semilogx plot.
% Shape determines the sharpness of the cut-off. 
% If shape =1, the filter is Gaussian; as shape increases 
% the filter shape becomes more and more rectangular.
% Declares global variables IFFDATA and IFFPARAMETERS
% Filtered signal is returned in the global variable IFFDATA(3,:,:)
% so declare IFFDATA as global before calling iff.m to
% have access to IFFDATA outside of iff.m
%  T. C. O'Haver (toh@umd.edu), version 3, June 2008
% Slider function by matthew jones (matt_jones@fastmail.fm), 2004

global IFFDATA % IFFDATA=[x;y;ry];
global IFFPARAMETERS % IFFPARAMETERS=[center width shape mode]
ry=ones(size(x));
IFFDATA=[x;y;ry];

close
signalstring='Original signal';

% Adjust x and y vector shape to 1 x n (rather than n x 1)
x=reshape(x,1,length(x));
y=reshape(y,1,length(y));

% Initial values of filter parameters
center=1;
width=1;
shape=2;
mode=0;  % mode=0 for band-pass filter, mode=1 for band-reject (notch) filter
IFFPARAMETERS=[center width shape mode];

% Plot the signal and its power spectrum
ry=RedrawFourierFilter(x,y,center,width,shape,mode);
h=figure(1);
h2=gca;
  
% Maximum ranges of the sliders (change as needed)
MaxCenter=(length(x)/2)-1;
MaxWidth=length(x);
MaxShape=20;

% Draw the sliders
rtslid(h,@iffcenter,h2,1,'Scale',[0 MaxCenter],'Def',center,'Back',[0.9 0.9 0.9],'Label','Center','Position',[0.03 0.24 0.03 0.68]);
rtslid(h,@iffmode,h2,0,'Scale',[0 3],'Def',.1,'Back',[0.9 0.9 0.9],'Label','Mode','Position',[0.03 0.05 0.03 0.1]);
rtslid(h,@iffwidth,h2,0,'Scale',[1 MaxWidth],'Def',width,'Back',[0.9 0.9 0.9],'Label','Width','Position',[0.95 0.24 0.03 0.68]);
rtslid(h,@iffshape,h2,0,'Scale',[1 MaxShape],'Def',shape,'Back',[0.9 0.9 0.9],'Label','Shape','Position',[0.95 0.05 0.03 0.1]);
ry=IFFDATA(3,:,:);

function iffcenter(n,h)
% Re-draws graph when the center slider is moved
% Tom O'Haver, May 2008
global IFFDATA % IFFDATA=[x;y;ry];
global IFFPARAMETERS % IFFPARAMETERS=[center width shape mode]
center=round(n);
IFFPARAMETERS(1)=center;
width=IFFPARAMETERS(2);
shape=IFFPARAMETERS(3);
mode=IFFPARAMETERS(4);
x=IFFDATA(1,:,:);
y=IFFDATA(2,:,:);
ry=RedrawFourierFilter(x,y,center,width,shape,mode);
IFFDATA(3,:,:)=ry;
axes(h);
h2=gca;

function iffmode(n,h)
% Re-draws graph when the width slider is moved
% Tom O'Haver, May 2008
global IFFDATA % IFFDATA=[x;y;ry];
global IFFPARAMETERS % IFFPARAMETERS=[center width shape mode]
mode=round(n);
center=IFFPARAMETERS(1);
width=IFFPARAMETERS(2);
shape=IFFPARAMETERS(3);
IFFPARAMETERS(4)=mode;
x=IFFDATA(1,:,:);
y=IFFDATA(2,:,:);
ry=RedrawFourierFilter(x,y,center,width,shape,mode);
IFFDATA(3,:,:)=ry;
axes(h);
h2=gca;  

function iffwidth(n,h)
% Re-draws graph when the width slider is moved
% Tom O'Haver, May 2008
global IFFDATA % IFFDATA=[x;y;ry];
global IFFPARAMETERS % IFFPARAMETERS=[center width shape mode]
width=round(n);
center=IFFPARAMETERS(1);
IFFPARAMETERS(2)=width;
shape=IFFPARAMETERS(3);
mode=IFFPARAMETERS(4);
x=IFFDATA(1,:,:);
y=IFFDATA(2,:,:);
ry=RedrawFourierFilter(x,y,center,width,shape,mode);
IFFDATA(3,:,:)=ry;
axes(h);
h2=gca;

function iffshape(n,h)
% Re-draws graph when the shape slider is moved
% Tom O'Haver, May 2008
global IFFDATA % IFFDATA=[x;y;ry];
global IFFPARAMETERS % IFFPARAMETERS=[center width shape mode]
shape=round(n);
center=IFFPARAMETERS(1);
width=IFFPARAMETERS(2);
IFFPARAMETERS(3)=shape;
mode=IFFPARAMETERS(4);
x=IFFDATA(1,:,:);
y=IFFDATA(2,:,:);
ry=RedrawFourierFilter(x,y,center,width,shape,mode);
IFFDATA(3,:,:)=ry;
axes(h);
h2=gca;

function ry=RedrawFourierFilter(xvector,yvector,centerfrequency,filterwidth,filtershape,mode)
% Separate graph windows for the original and filtered signals.
% Computes and plots fourier filter for signal yvector.  Centerfrequency
% and filterwidth are the center frequency and width of the
% pass band, in harmonics. 'filtershape' determines the sharpness of the 
% cut-off. If filtershape =1, the filter is Gaussian; as filtershape 
% increases the filter shape becomes more and more rectangular. 
% mode=0 for band-pass filter, linear plot; mode=1 for  band-pass filter,
% semilogx plot; mode=2 for band-reject (notch) filter, linear plot;
% mode=3 for band-reject (notch) filter, semilogx plot.
% In this version, the x-axis of the bottom window is expressed in Hz
% and CenterF and WidthF are the filter center and width in Hz.
%  T. C. O'Haver (toh@umd.edu),  version 2, May, 2008
global IFFDATA
global IFFPARAMETERS
fy=fft(yvector);
lft1=[1:(length(fy)/2)];
lft2=[(length(fy)/2+1):length(fy)];
% Compute filter shape
ffilter1=ngaussian(lft1,centerfrequency+1,filterwidth,filtershape);
ffilter2=ngaussian(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
ffilter=[ffilter1,ffilter2];
modestring='Band-pass:';
if mode>1, 
    ffilter=1-ffilter; 
    modestring='Band-reject (notch):';
end
if length(fy)>length(ffilter), ffilter=[ffilter ffilter(1)];,end
ffy=fy.*ffilter;  % Multiply filter by Fourier Transfork of signal
ry=real(ifft(ffy));

subplot(3,1,1)  % Plot original signal in top plot
plot(xvector,yvector,'b');
title(['Original Signal     x=sec.'])
axis([xvector(1) xvector(length(xvector)) min(yvector) max(yvector)]);  

subplot(3,1,2)  % Plot filtered signal in middle plot
plot(xvector,ry,'r');
title('Filtered signal    x=sec.')
axis([xvector(1) xvector(length(xvector)) min(ry) max(ry)]);  

subplot(3,1,3)    % Plot power spectrum and filter in lower plot
py=fy .* conj(fy); % Compute power spectrum
plotrange=1:length(fy)/2;
f=((plotrange-1)./range(xvector));
switch mode,
  case {0,2}
    plot(f,real(py(plotrange)),f,max(real(py)).*ffilter(plotrange),'r')
  case {1,3}
    semilogx(f,real(py(plotrange)),f,max(real(py)).*ffilter(plotrange),'r')
  otherwise,
end
title('BLUE = Power spectrum of signal      RED = Filter spectrum   x=Hz')
xlabel([ modestring '   CenterFreq = ' num2str(centerfrequency./range(xvector)) '    FreqWidth = ' num2str(filterwidth./range(xvector)) '     Shape =  ' num2str(filtershape)])
if centerfrequency+filterwidth/2<length(fy)/4,
    axis([0 max(f)/2 min(py) 1.1*max(real(py(plotrange)))]), 
else
    axis([0 max(f) min(py) 1.1*max(real(py(plotrange)))]), 
end


function g = ngaussian(x,pos,wid,n)
%  ngaussian(x,pos,wid) = peak centered on x=pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  Shape is Gaussian when n=1, becomes more rectangular as n increases.
% Example: ngaussian([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
g = exp(-((x-pos)./(0.6006.*wid)) .^(2*round(n)));

function h = rtslid(fig,f,hh,varargin)

%RTSLID Slider widget that responds to dragging realtime
%
%       h2 = RTSLID(h,@Fcn,h3) places a slider on the figure
%       with handle h, and returns the handle to the slider
%       in h2.  The second argument is a pointer to a
%       function that takes a two input arguments, which
%       should be the function you wish to be responding to
%       changes in the slider position.  The first of these
%       input arguments will be a number between 0 and 1
%       (if using the default scale) that corresponds to the
%       slider vertical position, while the second is the
%       handle of the axes to that responds to the slider
%       movement (as a result of Fnc). This handle should be
%       passed to RTSLID as the third input argument h3.
%
%       Alternatively instead of passing a function handle
%       as the second argument, a set of commands can be
%       passed as a string.  In order to use the slider
%       output value in such a set of commands, the lower
%       case variable 'out' should be used which contains
%       this number (Look at bottom of DEMO1.m for example).
%
%       h2 = RTSLID(h,@fcn,h3,style) allows the user to
%       specify the style of the slider by including an
%       optional third input argument.
%       The different styles affect the slider's response
%       to mouse movement, and are:
%
%               0 - (Default) Jumps to the mouse pointer
%                   on click, and locks to the mouse on drag.
%               1 - Only moves with mouse drags.
%
%       Additional parameters are passed as parameter pairs,
%       ie. rtslid(h,@fcn,h3,'param1',v1,'param2',v2,'...).
%       Inclusion of the fourth argument 'style' is optional.
%       These parameter names are listed below (all parameter
%       names are case insensitive) and all are optional:
%
%               'POSITION'  -   Control the position and size
%                               of the slider by passing a
%                               four element vector in the
%                               range [0 to 1] as a multiple
%                               of current figure size.
%                               Use [x y width height].
%                               {Default: [0.03 0.1 0.03 0.8]}
%               'BACK'      -   Either a string of the name
%                               of a valid colormap, or a 3-
%                               element RGB vector in the
%                               range [0 to 1].
%                               {Default: 'jet'}
%               'SCALE'     -   A two element vector of the
%                               lower and upper output limits
%                               for the slider.
%                               {Default: [0 1]}
%               'LABEL'     -   A string for the slider label.
%                               {Default: ''}
%               'DEF'       -   A default initial value (with
%                               respect to the SCALE), although
%                               does not cause output.
%                               {Default: 0.5}
%               'BUTMOT'    -   Additional mouse callback for
%                               WindowButtonMotionFcn of the
%                               slider figure, to allow other
%                               processes to use this callback.
%                               {Default: ''}
%               'MOUSEDOWN' -   Additional mouse callback for
%                               WindowButtonDownFcn of the
%                               slider figure, to allow other
%                               processes to use this callback.
%                               {Default: ''}
%               'MOUSEUP'   -   Additional mouse callback for
%                               WindowButtonDownFcn of the
%                               slider figure, to allow other
%                               processes to use this callback.
%                               {Default: ''}
%
%   matthew jones (matt_jones@fastmail.fm), 2004

DEFstyle = 0;
DEFback = 'jet';
DEFscale = [0 1];
DEFlabel = '';
DEFdef = 0.5;

style = DEFstyle;
back = DEFback;
scale = DEFscale;
label = DEFlabel;
def = DEFdef;

butmot = '';
mousedown = '';
mouseup = '';

if nargin>3,    % Set the style and any other parameters
    [style,args] = parseparams(varargin);
    if (length(style)==0),
        style = 0;
    else
        style = style{1};
    end
    try
        for l=1:(length(args)/2),
            inparam = args((l*2)-1);
            switch upper(inparam{1}),
            case 'POSITION' % Adjust the size and position of slider
                inparam = args(l*2);
                pos = inparam{1};
            case 'BACK'     % Change the background colour
                inparam = args(l*2);
                back = inparam{1};
            case 'SCALE'    % Change the output scale limits
                inparam = args(l*2);
                scale = inparam{1};
            case 'LABEL'    % Give the slider a label
                inparam = args(l*2);
                label = inparam{1};
            case 'DEF'      % Set the initial position
                inparam = args(l*2);
                def = inparam{1};
            case 'BUTMOT'   % Add extra WindowButtonMotionFcn commands
                inparam = args(l*2);
                butmot = inparam{1};
            case 'MOUSEDOWN'   % Add extra WindowButtonDownFcn commands
                inparam = args(l*2);
                mousedown = inparam{1};
            case 'MOUSEUP'   % Add extra WindowButtonUpFcn commands
                inparam = args(l*2);
                mouseup = inparam{1};
            end
        end
    catch
        error('Incorrect input arguments!');
    end
end

def2 = (def-scale(1))*(1000/(scale(2)-scale(1)));   % The initial position now in range [0 1000]

figure(fig);
params.fig = fig;   % Set so that focus returns to slider figure ofter plotting
params.butmot = butmot;
params.mousedown = mousedown;
params.mouseup = mouseup;

params2 = get(fig,'UserData');
if (length(params2)==0),    % If there are no sliders already,
                            % initialise the first elements of 'params'.
                            
    if exist('pos'),        % Draw the slider on its own axes
        h = axes('Position',pos);
    else
        h = axes('Position',[0.03 0.1 0.03 0.8]);
    end

	if ischar(back),        % Any string wll cause the background to be 'jet'
        try
            eval(['colormap(''',back,''');']);
        catch
            error('If passing a string for background, string must be valid colormap!');
        end
        back = repmat(linspace(0,1,1000).',[1,10]);
        imagesc(back);
    else                    % Or construct an RGB matrix m for the background
        m = ones(1000,10,3);
        for l=1:3,
            m(:,:,l) = back(l)*m(:,:,l);
        end
        image(m);
    end
            
	title(label);
    
    % Make the slider look nice
    set(h,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','','YDir','normal','LineWidth',2);
	set(fig,'DoubleBuffer','on');   % Double buffering required for smooth display

    % Draw the small black box that indicates current position
    ind = patch([0 0 10 10],[def2-20 def2+20 def2+20 def2-20],[0.25 0.25 0.25]);axis tight;

	params.mdown = 0;   % Initialise mousedown state
	params.ind = ind;   % Save the handle to the position indicator
	params.funhand{1} = f;  % Save the function handle or string
	params.h = h;       % Save the handle to this slider
	params.h2 = hh;     % Save the handle to the plotting axes (for this slider)
    if exist('scale'),  % Set the output limits
        params.scale{1} = scale;
    else
        params.scale{1} = DEFscale;
    end
    params.noslids = 1; % Keep tab of how many sliders in this figure
else
    if (params2.noslids<8), % If there is one or more sliders on this figure (and less
                            % than an upper limit), concatenate fields of 'params'
        params.noslids = params2.noslids+1; % Keep tab of how many sliders in this figure
        if exist('pos'),    % Draw the slider on its own axes
            h = axes('Position',pos);
        else
            h = axes('Position',[(((params.noslids-1)*0.06)+0.03) 0.1 0.03 0.8]);
        end
        
        if ischar(back),    % Any string wll cause the background to be 'jet'
            back = repmat(linspace(0,1,1000).',[1,10]);
            imagesc(back);
        else                % Or construct an RGB matrix m for the background
            m = ones(1000,10,3);
            for l=1:3,
                m(:,:,l) = back(l)*m(:,:,l);
            end
            image(m);
        end
    
        title(label);
        
        % Make the slider look nice
        set(h,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','','YDir','normal','LineWidth',2);
		set(fig,'DoubleBuffer','on');
	
        % Draw the small black box that indicates current position
        ind = patch([0 0 10 10],[def2-20 def2+20 def2+20 def2-20],[0.25 0.25 0.25]);axis tight;
 
        % Copy old parameters and concatenate new values
        params.mdown = [params2.mdown 0];   % Initialise new mousedown state
		params.ind = [params2.ind ind];     % Save the handle to the position indicator
		params.funhand = params2.funhand;   % Save the function handle or string
        params.funhand{params.noslids} = f;
		params.h = [params2.h h];           % Save the handle to this slider
		params.h2 = [params2.h2 hh];        % Save the handle to the plotting axes (for this slider)
        for l=1:(params.noslids-1),
            params.scale{l} = params2.scale{l};
        end
        if exist('scale'),  % Set the output limits
            params.scale{params.noslids} = scale;
        else
            params.scale{params.noslids} = DEFscale;
        end
        params.style = params2.style;       % Copy list of slider styles
        params.motfcn = params2.motfcn;     % Copy the WindowButtonMotionFcn list
        params.current = params2.current;   % Copy the 'current' slider positions
    end

end

if (style==0),
    params.style(params.noslids) = 0;               % Set the slider style
    params.motfcn{params.noslids} = {@butmotfcn0};  % Set the corresponding WindowButtonMotionFcn
else
    params.style(params.noslids) = 1;               % Set the slider style
    params.current(params.noslids) = def2;          % Initialise slider position (only required for style 1)
    params.motfcn{params.noslids} = {@butmotfcn1};  % Set the corresponding WindowButtonMotionFcn
end

% Now save these params and set the figure's mouse callback functions
set(fig,'UserData',params,'WindowButtonDownFcn',{@butdownfcn},'WindowButtonUpFcn',{@butupfcn});


% Callback functions defined here to keep any variables created
% localto these functions.
%------------------------------------------------------------------------%
function butdownfcn(in,varargin)

params=get(gcf,'UserData');     % Get current parameters
for l=1:params.noslids,
    p=get(params.h(l),'CurrentPoint');
    if ((p(1)>0) & (p(1)<10.7) & (p(3)>0) & (p(3)<1000)),   % Find out if mouse is on slider l
        params.mdown(l)=1;                                  % Set slider MouseDown sate to 1
        if (params.style(l)==0),                            % If jump-to-click style, use new
            params.old=p(3);                                % mouse position as output
            params.current(l)=p(3);                         % Save current position
            set(params.ind(l),'YData',[p(3)-20 p(3)+20 p(3)+20 p(3)-20]);   % Update indicator
            scal=params.scale{l};
            out=(p(3)*(scal(2)-scal(1))/1000);              % Convert output to desired scale
            out=out+scal(1);
        else                                                % Else use old position for output
            scal=params.scale{l};
            out=(params.current(l)*(scal(2)-scal(1))/1000); % Convert output to desired scale
            out=out+scal(1);
            params.old=p(3);
        end
        axes(params.h2(l));
        if ischar(params.funhand{l}),                       % If passed string as function
            eval(params.funhand{l});                        % evaluate using eval()
        else                                                % Else if passed function handle
            feval(params.funhand{l},out,params.h2(l));      % evaluate using feval()
        end
        figure(params.fig);
        set(gcf,'UserData',params,'WindowButtonMotionFcn',params.motfcn{l});    % Save new params and
    end;                                                                        % activate slider's Motion Function
end
if (length(params.mousedown)>0),    % If additional MouseDown parameters passed,
    eval(params.mousedown);         % evaluate these too.
end


%------------------------------------------------------------------------%
function butupfcn(in,varargin)

params=get(gcf,'UserData'); % Get old parameters
for l=1:params.noslids,
     params.mdown(l) = 0;   % Set ALL slider MouseDown states to zero
end
set(gcf,'UserData',params); % Save these new states

if (length(params.butmot)==0),  % If no additional ButMot commands passed,
    set(gcf,'WindowButtonMotionFcn','');    % Clear figure WindowButtonMotionFcn
else                            % Otherwise set it to these extra commands
    set(gcf,'WindowButtonMotionFcn',params.butmot);
end
    
if (length(params.mouseup)>0),  % Evaluate any additional MouseUp commands passed
    eval(params.mouseup);
end


%------------------------------------------------------------------------%
function butmotfcn0(in,varargin)
% Mouse motion function for style 0 (only active during mouse press)

params = get(gcf,'UserData'); % Get the saved parameters
for l=1:params.noslids,
   if (params.mdown(l)),
       p = get(params.h(l),'CurrentPoint');   % Get mouse position relative to the slider axis
       p = p(3);                              % Jump to mouse position
       p = max([0 min([1000 p])]);            % Constrain position
       if (abs(params.old-p)>0),            % Only act if mouse has moved
           scal = params.scale{l};            % Get the lower and upper limits
           out = (p*(scal(2)-scal(1))/1000);  % Convert from [0 to 1000] to new scale
           out = out+scal(1);
           axes(params.h2(l));              % Avoid plotting on the slider axis
           if ischar(params.funhand{l}),    % If string passed as function
               eval(params.funhand{l});     % Evaluate the string
           else                             % Else if function handle passed
               feval(params.funhand{l},out,params.h2(l));   % Perform function
           end
           set(params.ind(l),'YData',[p-20 p+20 p+20 p-20]);% Update position indicator
           figure(params.fig);              % Make sure focus returns to slider after plotting
       end
   end
end
params.old = p;             % Save old mouse position (vertical only)
set(gcf,'UserData',params); % Save all parameters
if (length(params.butmot)>0),   % Evaluate any extra ButMot functions passed to rtslid
    eval(params.butmot);
end


%------------------------------------------------------------------------%
function butmotfcn1(in,varargin)
% Mouse motion function for style 1 (only active during mouse press)

params = get(gcf,'UserData'); % Get the saved parameters
for l=1:params.noslids,
	if (params.mdown(l)),                       % If mouse was clicked on slider l
		p = get(params.h(l),'CurrentPoint');    % Get mouse position relative to the slider axis
		dy = p(3)-params.old;                   % Calculate change in mouse position
		if (abs(dy)>0),                         % If mouse has moved
			p2 = params.current(l)+dy;          % Find new slider position
			params.old = p(3);                  % Save this new mouse position
			p2=max([0 min([1000 p2])]);         % Constrain position
			params.current(l) = p2;             % Save current slider position
			scal = params.scale{l};             % Get the lower and upper limits
			out = (p2*(scal(2)-scal(1))/1000);  % Convert from [0 to 1000] to new scale
			out = out+scal(1);
			axes(params.h2(l));                 % Avoid plotting on the slider axis
			if ischar(params.funhand{l}),       % If string passed as function
			    eval(params.funhand{l});        % Evaluate the string
			else                                % Else if function handle passed
			    feval(params.funhand{l},out,params.h2(l));  % Perform function
			end
			set(params.ind(l),'YData',[p2-20 p2+20 p2+20 p2-20]);% Update position indicator
			figure(params.fig);                 % Make sure focus returns to slider after plotting
		end
	end
	set(gcf,'UserData',params); % Save all parameters
end
if (length(params.butmot)>0),   % Evaluate any extra ButMot functions passed to rtslid
    eval(params.butmot);
end