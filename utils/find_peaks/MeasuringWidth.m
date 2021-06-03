% Comparison of twho different meathod of measuring the width (full width
% at half maximum) of a peak. The peakfit.m and halfwidth.m function use
% two different methods of measuring width that ideally give exactly the
% same results. Ideally. Peakfit.m fits a Gaussian function to the data and
% reports the width of that Gaussian, and halfwidth uses an interpolation
% method to directly measure FWHM no matter what the shape of the peak. The
% two methods agree only if the data are a perfect Gaussian.
%
% The little script generates a perfect Gaussian (PureGaussian), then
% measures its width both ways.  There are three sources of "imperfection"
% that you can optionally add: 
% 
% (1) the Gaussian is sampled with an x-axis increment of "deltax", the
% smaller the better. 
% (2) random white noise with standard deviation "noise"  is added to the y
% data 
% (3) a flat baseline equal to "baseline" is added
%
% If, and only if, the deltax is very small, AND the noise and baseline are
% zero, AND if the peak is really Gaussian, then the two width measurement
% will be equal. Try changing these values in the first 4 lines. You can
% also change exp(-(x).^2) to exp(-(x).^4) to create a perfect
% NON-Gaussian; the widths disagree for that, too.
%
format compact
format short g
deltax=.001; % increment betwee x values
x=-5:deltax:5; % create x vector
noise=0.0001; % Optional random white noise 
baseline=0; % Optional baseline
PureGaussian=exp(-(x).^2); % Perfect Gaussian (change to exp(-(x).^4) for a perfect NON-Gaussian)
y=baseline+PureGaussian+noise.*randn(size(x)); % add baseline and noise
FWHM=halfwidth(x,y); % Measure FWHM by interpolation
P=peakfit([x;y]); % fit Gaussian to data, return peak table P
GW=P(:,4); % Extract peak width from the peak table P
disp(' ')
disp('       deltax       noise       baseline     GaussianW     FWHM')
disp([deltax noise baseline GW FWHM])
