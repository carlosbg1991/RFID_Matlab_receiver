function R2=trypolyplot(x,y)
% R2=trypolyplot(x,y) fits a series of polynomials to the data x,y, for
% polynomial orders 1 through length(x)-1, returns the R2 values in a
% vector, and plots polynimial orede vs R2 as a bar plot. Shows that, for
% any data, R2 -> 1 as order -> length(x)-1. Tom O'Haver, 2014
% Example:
% x = [1:12];y =cos(x);trypolyplot(x,y))
% Related functions: plotit.m, trydatatrans.m
warning off
R2=zeros(1,length(x)-1);
for k=1:length(x)-1,
    %sizepolyfit=size(polyfit(x,y,k))
    polycoeff=polyfit(x,y,k);
    R2(k)=RSquared(polycoeff,x,y);
end
bar(1:length(x)-1,R2)
xlabel('Polynomial Order')
ylabel('Coefficient of Determination (R2)')

function RS=RSquared(polycoeff,x,y)
    % Compute the correlation coefficient and R-Squared
    if IsOctave,
        cc=corr(polyval(polycoeff,x'),y');
        RS=cc.^2;
    else
        cc=corrcoef(polyval(polycoeff,x),y);
        RS=cc(2).^2;
    end %   if IsOctave,
    
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

