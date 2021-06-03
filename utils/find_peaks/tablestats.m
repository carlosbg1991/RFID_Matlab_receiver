function PS=tablestats(P,displayit)
% function PS=tablestats(P,displayit). Function to compute a table of
% summary statistics of the peak data in peak table 'P', such as generated
% by any of the functions that return a table of peaks with at least 4
% columns listing peak number, height, width, and area. Computes the peak
% intervals (the x-axis interval between adjacent detected peaks) and the
% maximum, minimum, average, and percent standard deviation of each, and
% optionally displaying the histograms of the peak intervals, heights,
% widths, and areas in figure window 2. The optional last argument
% "displayit" = 1 if the histograms are to be displayed, otherwize not.
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% T. C. O'Haver, 2015.  Version 1. Last revised October, 2015
%
% Example:
% x=[0:.1:1000];y=9+9.*cos(x)+randn(size(x));
% P=findpeaksplot(x,y,0,-1,11,19,3);
% tablestats(P,1);
%
% Related functions:
% peakstats.m, findvalleys.m, findpeaksL.m, findpeaksb.m, findpeaksplot.m,
% findpeaksG.m, findpeaksnr.m, findpeaksGSS.m, findpeaksLSS.m,
% findpeaksfit.m.
warning off
if nargin==1;displayit=0;end  % displayit=0 if not specified in 8th argument
sizeP=size(P);
if sizeP(1)>1,
    PP=P;
    %  "diff" is the vector of x-axis intervals between peaks.
    diff=zeros(size(1:sizeP(1)-1));
    for n=1:sizeP(1)-1;
        diff(n)=max(PP(n+1,2)-PP(n,2));
    end
    PS(1,:)=[max(diff);max(PP(:,3));max(PP(:,4));max(PP(:,5))];
    PS(2,:)=[min(diff);min(PP(:,3));min(PP(:,4));min(PP(:,5))];
    PS(3,:)=[mean(diff);mean(PP(:,3));mean(PP(:,4));mean(PP(:,5))];
    if IsOctave,
        PS(4,:)=[100.*stdev(diff)./mean(diff);100.*stdev(PP(:,3))./mean(PP(:,3));100.*stdev(PP(:,4))./mean(PP(:,4));100.*stdev(PP(:,5))./mean(PP(:,5))];
    else
        PS(4,:)=[100.*std(diff)./mean(diff);100.*std(PP(:,3))./mean(PP(:,3));100.*std(PP(:,4))./mean(PP(:,4));100.*std(PP(:,5))./mean(PP(:,5))];
    end
    disp('Peak Summary Statistics')
    disp( [ num2str(sizeP(1)) ' peaks listed in peak table' ] )
    disp('          Interval      Height      Width          Area')
    disp( [ 'Maximum    ' num2str(max(diff)) '       ' num2str(max(PP(:,3))) '       ' num2str(max(PP(:,4)))  '       ' num2str(max(PP(:,5)))  ])
    disp( [ 'Minimum    ' num2str(min(diff)) '       ' num2str(min(PP(:,3))) '       ' num2str(min(PP(:,4)))  '       ' num2str(min(PP(:,5)))  ])
    disp( [ 'Mean       ' num2str(mean(diff)) '       ' num2str(mean(PP(:,3))) '       ' num2str(mean(PP(:,4)))  '       ' num2str(mean(PP(:,5)))  ])
    if IsOctave,
        disp( [ '% STD      ' num2str(100.*stdev(diff)./mean(diff)) '       ' num2str(100.*stdev(PP(:,3))./mean(PP(:,3))) '       ' num2str(100.*stdev(PP(:,4))./mean(PP(:,4)))  '       ' num2str(100.*stdev(PP(:,5))./mean(PP(:,5)))  ])
    else
        disp( [ '% STD      ' num2str(100.*std(diff)./mean(diff)) '       ' num2str(100.*std(PP(:,3))./mean(PP(:,3))) '       ' num2str(100.*std(PP(:,4))./mean(PP(:,4)))  '       ' num2str(100.*std(PP(:,5))./mean(PP(:,5)))  ])
    end
    if displayit,
        figure(2);
        subplot(2,2,1);hist(diff);title('Histogram of intervals between peak positions')
        subplot(2,2,2);hist(PP(:,3));title('Histogram of peak heights')
        subplot(2,2,3);hist(PP(:,4));title('Histogram of Gaussian peak widths')
        subplot(2,2,4);hist(PP(:,5));title('Histogram of Gaussian peak areas')
    else
        PS=0;
    end
end
% ----------------------------------------------------------------------
function sd=stdev(a)
% Octave and Matlab compatible standard deviation function
sa=size(a);
reshape(a,1,length(a));
if sa(1)>sa(2),
  sd=std(a);
else
  sd=(std(a'));
end;

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
    % be cached in a separate variable.
    