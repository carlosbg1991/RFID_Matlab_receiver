function fs=FrequencySpectrum(DataMatrix,y)
% Returns the real part of the Fourier power spectrum of x,y as a matrix.
% Data may be separate x and y vectors or a single n by 2 matrix. Returns
% an n by 2 matrix with frequency and amplitude. To plot: plot(fs)
% Version 2, Octoiber 2019, Tom O'Haver toh@umd.edu
%
% Example 1: 
% Frequency spectrum of noisy sine wave; data as separate x and y vectors
% x=0:.01:2*pi;
% y=sin(200*x)+randn(size(x));
% fs=FrequencySpectrum(x,y);
% plot(fs)
%
% Example 2: 
% As above, with data in single n by 2 matrix as input argument.
% x=0:.01:2*pi;
% y=sin(200*x)+randn(size(x));
% DataMatrix=[x;y]; %  314x2 matrix
% fs=FrequencySpectrum(DataMatrix); 
% plot(fs)
%
switch nargin % Process arguments
    % 'nargin' is the number of arguments
    case 1  % One argument only
        % Might be isignal(DataMatrix) ot isignal(Y-vector)
        % If data is in the wrong transposition, fix it.
        datasize=size(DataMatrix);  
        if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
        datasize=size(DataMatrix);
        if datasize(2)==1 %  Must be isignal(Y-vector)
            x=1:length(DataMatrix); % Create an independent variable vector
            y=DataMatrix;
        else
            % Must be isignal(DataMatrix)
            x=DataMatrix(:,1); % Split matrix argument
            y=DataMatrix(:,2);
        end
    case 2
        % Two arguments, might be separate x and y data vectors,
        % or one data matrix and a peak density estimate.
        if isscalar(y) % if second argument is scalar
            % Must be isignal(DataMatrix,y)
            % If DataMatrix is in the wrong transposition, fix it.
            datasize=size(DataMatrix);
            x=DataMatrix(:,1); % Split matrix argument
            y=DataMatrix(:,2);
        else % if second argument is not scalar
            % Must be isignal(x,y)
            xdatasize=size(DataMatrix);
            if xdatasize(1)<xdatasize(2),DataMatrix=DataMatrix';end
            x=DataMatrix;  % First argument is X
        end  % if isscalar
end
% 
x=reshape(x,[1 length(x)]);
y=reshape(y,[1 length(x)]);
fy=fft(y);
sy=fy .* conj(fy);
plotrange=1:length(fy)/2;
realsy=real(sy(plotrange));
f=((plotrange-1)./range(x));
fs=([f;realsy])';

function r=range(x)
r=max(x)-min(x);
