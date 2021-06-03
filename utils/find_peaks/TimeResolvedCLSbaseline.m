% Time Resolved Classical Least Squares (TRCLS) with baseline. Simulation
% of the separation of three components plus baseline by chromatograpy with
% an array detector. Change the unknown mixture concentrations in line 29,
% the random detector noise in line 18, the spectra in lines 33-48, and the
% chromatographic peaks in lines 62-84. Initial values (lines 17-57) are
% based on the three acetophenone positional isomers: o-methyl (o-MAP),
% m-methyl (m-MAP), and p-methyl (p-MAP), estimated from
% https://solutions.shimadzu.co.jp/an/n/en/hplc/jpl217011.pdf Note. This
% script requires the following downloadable functions in the Matlab path,
% all freely available at http://tinyurl.com/cey8rwh modelpeaks.m,
% expgaussian.m, measurepeaks.m, halfwidth.m, val2ind.m, plotit.m
clear
format compact
format short g

% -------Initial settings, to change as you wish (lines 18-55) --------
wincrement=1.5; % wavelength increment (= slit width in Shimadzu specs)
w=210:wincrement:320; % wavelength axis 
t=4.7:.0035:6; % Time axis (1.5 minutes total duration)
tau=10; % Exponential time constant (larger = more asymmetrical)
DetNoise=.002; % Detector noise standard deviation (Shimadzu specification)
SmoothWidth=5; % Smoothing of chromatographic peaks to reduce noise
%
% SOLUTION concentrations
% Known concentrations of components in the standard solution [1 3 2]
M=[400 400 400 0]; % The 4th number is the baseline amplitude
% Unknown concentrations of components in the simulated unknown [1 3 2] 
U=[.01 .01 .01 0.1]; % NOTE: The 4th number is the baseline amplitude
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
SpecBaseline=U(4).*ones(size(w)); % Spectrum of basline is flat
%
% CHROMATOGRAPHY of Component 1
ort=5.0; % ortho retention time (min) 
ow=.11 ; % ortho peak width (min)
oh=220 ; % ortho peak height (mAU)
% Chromatograpy of Component 2
mrt=5.2; % meta retention time (min) 
mw=.12 ; % meta peak width (min)
mh=220 ; % meta peak height (mAU)
% Chromatograpy of Component 3
prt=5.4; % para retention time (min)
pw=.13 ; % para peak width (min)
ph=260 ; % para peak height (mAU)
% Chromatographic baseline is modeled as a sloped straight line
baseline=((t-min(t)).*ones(size(t)));

% ------------ Calculations ---------------------------
% Each spectrum is modeled as the sum of 3 Gaussians + random detector noise
% Compute spectrum of component 1 (OrthoMAP)
oMAPspect=modelpeaks(w,3,1,omAU,opeaks,owidth,0,0)+DetNoise.*randn(1,length(w));
% Compute spectrum of component 2 (MetaMAP)
mMAPspect=modelpeaks(w,3,1,mmAU,mpeaks,mwidth,0,0)+DetNoise.*randn(1,length(w));
% Compute spectrum of component 3 (ParaMAP)
pMAPspect=modelpeaks(w,3,1,pmAU,ppeaks,pwidth,0,0)+DetNoise.*randn(1,length(w));
MixtureSpectrum=U(1)*oMAPspect+U(2)*mMAPspect+U(3)*pMAPspect+DetNoise.*randn(1,length(w));

figure(1)
clf
% Plot all three spectra in differenet colors
plot(w,oMAPspect,'b',w,mMAPspect,'g',w,pMAPspect,'r',w,SpecBaseline,'c')
xlabel('Wavelength, nm')
ylabel('mAU')
title ('Known spectra of 1 (blue), 2 (green), 3 (red), and baseline (cyan), matrix A')
axis([210 320 -20 410])
grid

% ---Chromatography----
% Chromatographic peaks (modeled as Exponentially modified Gaussians)
% Peak of component 1
oMAPchrom=oh*expgaussian(t,ort,ow,-tau); % tau is exponential time constant
oMAParea=trapz(t,oMAPchrom); % area by trapezoidal numerical integration
hwo=halfwidth(t,oMAPchrom);
rto=t(val2ind(oMAPchrom,max(oMAPchrom))); % Retention time
% Peak of component 2
mMAPchrom=mh*expgaussian(t,mrt,mw,-tau);% tau is exponential time constant
mMAParea=trapz(t,mMAPchrom); % area by trapezoidal numerical integration
hwm=halfwidth(t,mMAPchrom);
rtm=t(val2ind(mMAPchrom,max(mMAPchrom))); % Retention time
% Peak of component 3
pMAPchrom=ph*expgaussian(t,prt,pw,-tau);% tau is exponential time constant
pMAParea=trapz(t,pMAPchrom); % area by trapezoidal numerical integration
hwp=halfwidth(t,pMAPchrom);
rtp=t(val2ind(pMAPchrom,max(pMAPchrom))); % Retention time
BGarea=trapz(t,baseline); % Area under baseline
% Note: measured baseline will be subtrtacted in line 219

TrueRetentionTimes=[rto rtm rtp];
% Plot three chromatographic peaks in different colors
figure(2)
clf
chrom=oMAPchrom+mMAPchrom+pMAPchrom; % Sum of the three peaks
plot(t,oMAPchrom,'b',t,mMAPchrom,'g',t,pMAPchrom,'r',t,chrom,'k.')
ylabel('mAU')
xlabel('Time, minutes')
title('Chromatographic peaks of separate component standard, matrix H')
maxy=1.05*max(chrom); 
text(ort,maxy,'1');text(mrt,maxy,'2');text(prt,maxy,'3')
% Add baseline
BG=ones(size(oMAPspect));
A=[oMAPspect;mMAPspect;pMAPspect;BG]+DetNoise.*randn(1,length(w));
H=[oMAPchrom';mMAPchrom';pMAPchrom';baseline]+DetNoise.*randn(1,length(t));
E=A'./M; % Analytical sensitivity

% Alternative calculation #1, lines 113-146, to animate the evolution
%  of the chromatographic profiles in the subplots of Figure 3.
% C=zeros(3,length(t))';
% for x=1:length(t)
%     Do=(oMAPchrom(x)*oMAPspect)./(M(1)); % oMAP at time x
%     Dm=(mMAPchrom(x)*mMAPspect)./(M(2)); % mMAP at time x
%     Dp=(pMAPchrom(x)*pMAPspect)./(M(3)); % pMAP at time x
%     % Use unknown concentrations U to weight the
%     % contribution of each component in the spectrum
%     % and add DetNoise mAU random noise.
%     % Spectral profile at the detector at time t:
%     D=U(1)*Do + U(2)*Dm + U(3)*Dp + DetNoise.*randn(1,length(w)); 
%     % Use analytical sensitivity E to compute
%     % concentrations C of each component at time t(x)
%     C(x,:)=D/E'; % C = concentrations of each component at time t(x)
%     % Plot the chromatographic profile of each separate component
%     figure(3)
%     subplot(2,2,1)
%     y1=fastsmooth(C(:,1),SmoothWidth,3,1);
%     plot(t,y1,'b')
%     ylabel('mAU')
%     xlabel('Time, minutes')
%     title('Component 1')
%     subplot(2,2,2)
%     y2=fastsmooth(C(:,2),SmoothWidth,3,1);
%     plot(t,y2,'g')
%     ylabel('mAU')
%     xlabel('Time, minutes')
%     title('Component 2')
%     subplot(2,2,3)
%     y3=fastsmooth(C(:,3),SmoothWidth,3,1);
%     plot(t,y3,'r')
%     ylabel('mAU')
%     xlabel('Time, minutes')
%     title('Component 3')
% end

% Alternative calculation #2, lines 150-170, matrix math method
% rather than a loop: shorter and faster, but without animation
C=(H'.*U)+DetNoise.*randn(1,length(t))';
% Plot the chromatographic profile of each separate component
figure(3)
subplot(2,2,1)
y1=fastsmooth(C(:,1),SmoothWidth,3,1);
plot(t,y1,'b')
ylabel('mAU')
xlabel('Time, minutes')
title('Component 1')
subplot(2,2,2)
y2=fastsmooth(C(:,2),SmoothWidth,3,1);
plot(t,y2,'g')
ylabel('mAU')
xlabel('Time, minutes')
title('Component 2')
subplot(2,2,3)
y3=fastsmooth(C(:,3),SmoothWidth,3,1);
plot(t,y3,'r')
ylabel('mAU')
xlabel('Time, minutes')
title('Component 3')
subplot(2,2,4)
y4=fastsmooth(C(:,4),SmoothWidth,3,1);
plot(t,y4,'c')  % baseline in cyan
ylabel('mAU')
xlabel('Time, minutes')
title('Baseline')

CS=C*A; % The Chromatographic-Spectroscopic matrix, CS
% The CS matrix contains both the chromatographic (C) and Spectrocopic (S)
% data in one matrix. To dispay a 3-D surface or a contour plot, type
% "mesh(CS)" or "contour(CS)".
% To plot the spectrum at time "time" (in minutes), type 
% "plot(w,CS(val2ind(t,time),:))".
% To plot the chromatogram at any single wavelength "nm", type
% "plot(t,CS(:,val2ind(w,nm)))".

% Report and plot results
ComponentsDetected=rank(H'*A)
HalfWidths=[halfwidth(t,y1) halfwidth(t,y2) halfwidth(t,y3)];
to1=t(val2ind(y1,max(y1))); % Retention time 1 measured by CLS
to2=t(val2ind(y2,max(y2))); % Retention time 2 measured by CLS
to3=t(val2ind(y3,max(y3))); % Retention time 3 measured by CLS

disp('     Component 1  Component 2  Component 3   Baseline')
TrueRetentionTimes
CLSRetentionTimes=[to1 to2 to3]
disp(' ')
figure(4)
yobsd=[oMAPchrom mMAPchrom pMAPchrom baseline']*U'+DetNoise.*randn(1,length(t))';
yobsd=fastsmooth(yobsd,SmoothWidth,3,1);
plot(t,yobsd,'k',t,y1,'b',t,y2,'g',t,y3,'r',t,y4,'c')
ylabel('mAU')
xlabel('Time, minutes')
title('Chromatogram of mixture, matrix C')
maxy=.9*max(y1+y2+y3);
text(ort,maxy,'1');text(mrt,maxy,'2');text(prt,maxy,'3')
TruePeakAreas=[oMAParea mMAParea pMAParea BGarea].*U
CLSPeakAreas=[trapz(t,y1) trapz(t,y2) trapz(t,y3) trapz(t,y4)]
CLSPercentError=100.*[(CLSPeakAreas(1)-TruePeakAreas(1))./CLSPeakAreas(1) (CLSPeakAreas(2)-TruePeakAreas(2))./CLSPeakAreas(2) (CLSPeakAreas(3)-TruePeakAreas(3))./CLSPeakAreas(3)]
disp(' ')
% disp('Attempt to detect and measure peaks directly using measurepeaks.m')
try
    P=measurepeaks(t',yobsd-y4,0.00000001,0.1,11,11,0); % Subtract baseline
    sizeP=size(P);
    if sizeP(1)==3
        disp('Perpendicular Drop Method, baseline subtracted, using measurepeaks.m')
        PerpDropRetentionTimes=[P(1,2) P(2,2) P(3,2)]
        PerpDropAreas=[P(1,5) P(2,5) P(3,5)]
        PerpDropPercentError=-100.*[(TruePeakAreas(1)-PerpDropAreas(1))./TruePeakAreas(1) (TruePeakAreas(2)-PerpDropAreas(2))./TruePeakAreas(2) (TruePeakAreas(3)-PerpDropAreas(3))./TruePeakAreas(3)] 
        MeanAbsCLSerror=mean(abs(PercentError));
        MeanAbsPDerror= mean(abs(PerpDropPercentError));
        disp(['Area measurement by CLS better that PD by factor of ' num2str(MeanAbsPDerror/MeanAbsCLSerror)])
    end
    if sizeP(1)==2
        disp('Only 2 peaks detected by measurepeaks.m. ')
    end
    if sizeP(1)==1
        disp('Only 1 peak detected by measurepeaks.m.')
    end   
    if sizeP(1)>3
        disp('More than 3 peaks detected by measurepeaks.m.')
    end  
catch
    % disp('Perpendicular drop method failed.')
end


%  Un-comment the next two lines to see the output of measurepeaks.m
%     disp('          Peak      Position    PeakMax    Peak-valley    Perp drop   Tan skim')
%     disp(P)

% Un-comment the following lines to see a 3D mesh representation of the C*A
% matrix in figure 5.
% figure(5)
% mesh(C*A) 
% xlabel('Wavelength') 
% ylabel('Time, minutes')
% zlabel('mAU')