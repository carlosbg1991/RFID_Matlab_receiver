% Script that compares statistics of TFit method to single-wavelength,% simple regression, and weighted regression methods. Repeats any number of% times, computes mean, standard deviation, and deviation from true% absorbance (which you can set in line 22).% You can change the parameters in lines 22-29.% Required functions: gaussian, lorentzian, fitM% Example: Just type TFitStats in the command window.% Tom O'Haver, August 27, 2006. Modified for Octave, Oct. 2012clearformat compactformat short gglobal z c% Arrays:% y = true transmission spectrum, without noise or broadening% InstFunction = instrument function (modeled as a Gaussian)% yobsd = noisy instrumentally broadened spectrum% f = frequency coordinate vector% x = wavelength coordinate vector% Initial values of the user-adjustable parameters:absorbance=1  % True absorbance of analyte at line center (try 0.001 to 1000)repeats=50     %  Number of repeat calculations with new noise sampleArrayLength=128   % Number of points in signalwidth=10;     % FWHM of absorption peak InstWidth=20  % FWHM of broadening function (spectrometer slit width)noise=0.01    % Random noise level when InstWidth = 1straylight=.01 % May be a scalar or a vector of length ArrayLength (Slider adjustable)IzeroShift=.01 % Random shifts in the 100% T intensity due to background absorption% Define frequency and wavelength coordinate vectorsx=[1:ArrayLength]';j=[-ArrayLength/2:(ArrayLength/2)-1]';f=(ArrayLength/2)-abs(j);% Calculate noisy instrumentally-broadened transmission profile, yobsd,% by convoluting true transmission spectrum y with instrument function% InstFunction and adding noise.  % Note:  To model gaussian absorption, change 'lorentzian' to 'gaussian'TrueSpectrum=lorentzian(x,(ArrayLength/2),width); y=10 .^ (absorbance .* (-TrueSpectrum));fy=fft(y);InstFunction=gaussian(f,0,InstWidth);  % define Gaussian instrument function centered on zerofa=fft(InstFunction);fy1=fy.*fa;            % Convolve the transmisison profile with the instrument function (InstWidth) yobsd=real(ifft(fy1));  % by multiplying Fourier transforms and inverse transforming the result.yo=yobsd./sum(InstFunction);for k=1:repeats, % Repeat k times with different random noise samples  yobsd=straylight+yo+((noise/InstWidth).*randn(size(yo))).*sqrt(yo);   % Add simulated photon noise  yobsd=yobsd.*(1-straylight);  yobsd=yobsd.*(1+IzeroShift.*randn); % Random shifts in Izero from sample to sample after instrument is zeroed  % Conventional methods  SingleWavelengthAbsorbance=-log10(yobsd(ArrayLength./2));  SimpleRegression=TrueSpectrum\(-log10(yobsd));  Background=ones(size(y));  weight=y;  WeightedRegression=([weight weight] .* [Background TrueSpectrum])\(-log10(yobsd) .* weight);  % Curve fitting methodstart=10; % Because of the very large dynamic range of absorbance, two start values are if SingleWavelengthAbsorbance<1,start=1;end  % used to prevent stalling on local optima.  lam=fminsearch(@(lambda)(fitM(lambda,yobsd,TrueSpectrum,InstFunction,straylight)),start);    results(k,:)=[absorbance SingleWavelengthAbsorbance SimpleRegression WeightedRegression(2) lam];enddisp('          True A    SingleW    SimpleR      WeightR      TFit')  MeanResult=mean(results)PercentRelativeStandardDeviation=100.*(std(results)./abs(MeanResult))PercentAccuracy=100.*(MeanResult-absorbance)./absorbancedisp('PercentAccuracy is the percent deviation of the mean from the true absorbance')% (Optional) Plot spectral profiles (uncomment next 6 lines)% plot(x,real(yobsd),'r.',x,real(y),'g',x,real(c)*z,'b',x,gaussian(x,ArrayLength/2,InstWidth),'m:'); % text(5,1.32,'Green = Reference spectrum      Dotted Magenta = Instrument function'); % text(5,1.25,'                    Red = Observed T     Blue = Fit to observed T'); % xlabel('Wavelength'); ylabel('Transmission');% title(['True absorbance = ' num2str(absorbance) '    Abs.Width = ' num2str(round(10*width)/10)  '   Inst.Width = ' num2str(round(10*InstWidth)/10) '    straylight= ' num2str(round(1000*mean(straylight))/10) '%' ]);% axis([0 ArrayLength 0 1.4]);% function err = fitM(lam,yobsd,Spectra,InstFun,StrayLight)% % Fitting function for broadened absorption of any number of components% % yobsd =  observed transmission spectrum (column vector)% % Sprecta = reference spectra for each component, one component/column% % InstFunction = Instrument function or slit function. (column vector)% % StrayLight = fractional stray light (scalar or column vector)% % Typical use: FMINSEARCH('fitM',start,options,yobsd,Spectra,InstFunction,StrayLight)% % yobsd, Spectra, and InstFunction must have same number of rows (wavelengths)% %  T. C. O'Haver, August 2006% global z% global c% A = StrayLight + (10 .^ -(Spectra*lam'));% fy=fft(A);% fa=fft(InstFun);% fy1=fy.*fa;                % z=real(ifft(fy1))./sum(InstFun);   % c = z\yobsd;% q = z*c;% err = norm(q-yobsd);