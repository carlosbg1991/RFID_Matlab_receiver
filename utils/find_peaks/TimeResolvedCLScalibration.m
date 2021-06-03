% Calibaration curves of Time Resolved Classical Least Squares (TRCLS).
% Simulation of the separation of three components by chromatograpy with an
% array detector. You can change the concentrations of the unknown mixture
% in line 29, the random detector noise in line 18, the spectra in lines
% 33-48, and the chromatographic peaks in lines 62-84. To make it more
% realistic, simulated sample spectra and chromatographic peaks are based
% on the three acetophenone positional isomers: o-methyl (o-MAP), m-methyl
% (m-MAP), and p-methyl (p-MAP), estimated from the graphs in the technical
% report from Shimadazu:
% https://solutions.shimadzu.co.jp/an/n/en/hplc/jpl217011.pdf Note. This
% script requires the following downloadable functions in the Matlab path,
% all freely available at http://tinyurl.com/cey8rwh modelpeaks.m,
% expgaussian.m, measurepeaks.m, plotit.m, val2ind.m
clear
format compact
format short g

% -------Initial settings, to change as you wish (lines 19-56) --------
NumTrials=10; % Number of different concentrations on the calibration curves
MaxConc=.001; % Maximum concentration on the calibration curves
wincrement=1.5; % wavelength increment (= slit width in Shimadzu specs)
w=210:wincrement:320; % wavelength axis
t=4.7:.0035:6.2; % Time axis (1.5 minutes total duration)
tau=21; % Exponential time constant (larger = more asymmetrical)
DetNoise=.002; % Detector noise standard deviation (Shimadzu specification)
SmoothWidth=21; % Smoothing of chromatographic peaks to reduce noise
%
% SOLUTION concentrations (Note: unknowns are set in line 58)
% Known concentrations of components in the standard solution [1 3 2]
M=[400 400 400];
%
% SPECTRAL characteristics of component 1 (based on OrthoMAP)
opeaks=[200 245 290]; % Peak wavelength, nm
omAU=[900 220 40];  % Peak absorbance, mAU
owidth=[20 25 30];  % Peak width, nm (FWHM)
% Spectral characteristics of component 2 (based on MetaMAP)
mpeaks=[200 250 290]; % Peak wavelength, nm
mmAU=[900 260 40]; % Peak absorbance, mAU
mwidth=[20 25 35];% Peak width, nm (FWHM)
% Spectral characteristics of component 3 (based on ParaMAP)
ppeaks=[200 255 290]; % Peak wavelength, nm
pmAU=[500 360 4]; % Peak absorbance, mAU
pwidth=[18 30 60];% Peak width, nm (FWHM)
%
% CHROMATOGRAPHY of Component 1
ort=5.0; % ortho retention time (min) 
ow=.11 ; % ortho peak width (min)
oh=220 ; % ortho peak height (mAU)
% Chromatograpy of Component 2
mrt=5.3; % meta retention time (min) 
mw=.12 ; % meta peak width (min)
mh=220 ; % meta peak height (mAU)
% Chromatograpy of Component 3
prt=5.6; % para retention time (min)
pw=.13 ; % para peak width (min)
ph=260 ; % para peak height (mAU)


% ------------ Calculations ---------------------------
for trial=1:NumTrials
    % Unknown concentrations of components in the simulated unknown [1 3 2]
    U=[MaxConc*rand MaxConc*rand MaxConc*rand]+MaxConc/32;
    ConcMatrix(trial,:)=U;
    % Each spectrum is modeled as the sum of 3 Gaussians + random detector noise
    % Compute spectrum of component 1 (OrthoMAP)
    oMAPspect=modelpeaks(w,3,1,omAU,opeaks,owidth,0,0)+DetNoise.*randn(1,length(w));
    % Compute spectrum of component 2 (MetaMAP)
    mMAPspect=modelpeaks(w,3,1,mmAU,mpeaks,mwidth,0,0)+DetNoise.*randn(1,length(w));
    % Compute spectrum of component 3 (ParaMAP)
    pMAPspect=modelpeaks(w,3,1,pmAU,ppeaks,pwidth,0,0)+DetNoise.*randn(1,length(w));
    MixtureSpectrum=U(1)*oMAPspect+U(2)*mMAPspect+U(3)*pMAPspect+DetNoise.*randn(1,length(w));
   
    % ---Chromatography----
    % Chromatographic peaks (modeled as Exponentially modified Gaussians)
    % Peak of component 1
    oMAPchrom=oh*expgaussian(t,ort,ow,-tau); % tau is exponential time constant
    oMAParea=trapz(t,oMAPchrom); % area by trapezoidal numerical integration
    % hwo=halfwidth(t,oMAPchrom);
    rto=t(val2ind(oMAPchrom,max(oMAPchrom))); % Retention time
    % Peak of component 2
    mMAPchrom=mh*expgaussian(t,mrt,mw,-tau);% tau is exponential time constant
    mMAParea=trapz(t,mMAPchrom); % area by trapezoidal numerical integration
    % hwm=halfwidth(t,mMAPchrom);
    rtm=t(val2ind(mMAPchrom,max(mMAPchrom))); % Retention time
    % Peak of component 3
    pMAPchrom=ph*expgaussian(t,prt,pw,-tau);% tau is exponential time constant
    pMAParea=trapz(t,pMAPchrom); % area by trapezoidal numerical integration
    % hwp=halfwidth(t,pMAPchrom);
    rtp=t(val2ind(pMAPchrom,max(pMAPchrom))); % Retention time
    
    TrueRetentionTimes=[rto rtm rtp];
    chrom=oMAPchrom+mMAPchrom+pMAPchrom; % Sum of the three peaks
    A=[oMAPspect;mMAPspect;pMAPspect]+DetNoise.*randn(1,length(w));
    H=[oMAPchrom';mMAPchrom';pMAPchrom']+DetNoise.*randn(1,length(t));
    E=A'./M; % Analytical sensitivity
    % Alternative calculation #2, lines 150-170, matrix math method
    % rather than a loop: shorter and faster, but without animation
    C=(H'.*U)+DetNoise.*randn(1,length(t))';
    
    % Report and plot results
    ComponentsDetected=rank(H'*A);
    yobsd=[oMAPchrom mMAPchrom pMAPchrom]*U'+DetNoise.*randn(1,length(t))';
    yobsd=fastsmooth(yobsd,SmoothWidth,3);
    y1=fastsmooth(C(:,1),SmoothWidth,3);
    y2=fastsmooth(C(:,2),SmoothWidth,3);
    y3=fastsmooth(C(:,3),SmoothWidth,3);
    TruePeakAreas=[oMAParea mMAParea pMAParea].*U;
    CLSPeakAreas=[trapz(t,y1) trapz(t,y2) trapz(t,y3)];
    AreaMatrix(trial,:)=CLSPeakAreas;
    CLSPercentError(trial,:)=100.*[(CLSPeakAreas(1)-TruePeakAreas(1))./CLSPeakAreas(1) (CLSPeakAreas(2)-TruePeakAreas(2))./CLSPeakAreas(2) (CLSPeakAreas(3)-TruePeakAreas(3))./CLSPeakAreas(3)];
    % disp(' ')
    % disp('Attempt to detect and measure peaks directly using measurepeaks.m')
    try
        P=measurepeaks(t',yobsd,0.00000001,0.001,11,11,0);
        sizeP=size(P);
        if sizeP(1)==3
            % disp('Perpendicular Drop Method, using measurepeaks.m')
            PerpDropRetentionTimes=[P(1,2) P(2,2) P(3,2)];
            PerpDropAreas=[P(1,5) P(2,5) P(3,5)];
            PerpDropAreaMatrix(trial,:)=PerpDropAreas;
          
            PerpDropPercentError(trial,:)=-100.*[(TruePeakAreas(1)-PerpDropAreas(1))./TruePeakAreas(1) (TruePeakAreas(2)-PerpDropAreas(2))./TruePeakAreas(2) (TruePeakAreas(3)-PerpDropAreas(3))./TruePeakAreas(3)];
            MeanAbsCLSerror=mean(abs(PercentError));
            MeanAbsPDerror= mean(abs(PerpDropPercentError));
            %       disp(['Area measurement by CLS better that PD by factor of ' num2str(MeanAbsPDerror/MeanAbsCLSerror)])
        end
        if sizeP(1)==2
            %        disp('Only 2 peaks detected by measurepeaks.m. ')
        end
        if sizeP(1)==1
            %      disp('Only 1 peak detected by measurepeaks.m.')
        end
        if sizeP(1)>3
            %      disp('More than 3 peaks detected by measurepeaks.m.')
        end
    catch
        % disp('Perpendicular drop method failed.')
    end
    
end

% PLot calibration curves
figure(1)
clf
plotit(ConcMatrix(:,1),PerpDropAreaMatrix(:,1),1,'o');
title('Perpendicular drop, Component 1')
xlabel('Concentration,  ug/mL')
ylabel('Peak area, mAU-min')
figure(2)
clf
plotit(ConcMatrix(:,2),PerpDropAreaMatrix(:,2),1,'o');
title('Perpendicular drop Component 2')
xlabel('Concentration,  ug/mL')
ylabel('Peak area, mAU-min')
figure(3)
clf
plotit(ConcMatrix(:,3),PerpDropAreaMatrix(:,3),1,'o');
title('Perpendicular drop Component 3')
xlabel('Concentration,  ug/mL')
ylabel('Peak area, mAU-min')
PerpDropAverageAbsolutePercentError=mean(abs(PerpDropPercentError))

figure(4)
clf
plotit(ConcMatrix(:,1),AreaMatrix(:,1),1,'o');
title('CLS Component 1')
xlabel('Concentration,  ug/mL')
ylabel('Peak area, mAU-min')
figure(5)
clf
plotit(ConcMatrix(:,2),AreaMatrix(:,2),1,'o');
title('CLS Component 2')
xlabel('Concentration,  ug/mL')
ylabel('Peak area, mAU-min')
figure(6)
clf
plotit(ConcMatrix(:,3),AreaMatrix(:,3),1,'o');
title('CLS Component 3')
xlabel('Concentration,  ug/mL')
ylabel('Peak area, mAU-min')
CLSAverageAbsolutePercentError=mean(abs(CLSPercentError))
