% Asymmetrical Area Measurement Test 2. Compares standard deviations of
% peak area measurements for an asymmetrical peak,comparing (A) Gaussian
% estimation,(B) triangulation, (C) perpendicular drop method, and curve
% fitting by (D) exponentially braodened Gaussian, and (E) two overlapping
% Gaussians. Must have the following functions in the Matlab/Octave path:
% gaussian.m, expgaussian.m, findpeaksG.m, findpeaksT.m,
% autopeaks.m,and peakfit.m.
clf
disp(' ')
noise=0.02; % Increase noise to test for noise robustness
NumTrials=50; % Increase NumTrials to improve precision of standard deviations

xx=1:100;
for trial=1:NumTrials,
    % Uncomment one of these lines,15,17,18,or 21, to activate that shape.
    % Alternative shape: Two closely overlapping Guassians
    yy=2.*gaussian(xx,max(xx)/2.4,max(xx)/4)+gaussian(xx,max(xx)/2,max(xx)/3);
    % Alternative shape: exponentially-broadened Gaussian
    % yy=2.*expgaussian(xx,max(xx)/2.4,max(xx)/4,-7)';
    % Alternative shape: triangle
    % yy=2.*triangle(xx,max(xx)/2.4,max(xx)/4);
    % Alternative shape: single Gaussian
    % yy=2.*gaussian(xx,max(xx)/2.4,max(xx)/4);
    
    yy=yy+noise.*randn(size(xx));
    TotalArea(trial)=sum(yy)./(xx(2)-xx(1));
    
    % Gaussian estimation method (findpeaks.m)
    Pge=findpeaksG(xx,yy,.01,1,3,max(xx)/4,1);
    AreaGge(trial)=Pge(5)./(xx(2)-xx(1));
    hold on 
    ge=Pge(3).*gaussian(xx,Pge(2),Pge(4));
    hold off
    
    % Triangulation method (shown in upper-right panel)
    Ptm=findpeaksT(xx,yy,.01,1,3,5,1);
    AreaTm(trial)=Ptm(5)./(xx(2)-xx(1));
    
    % Perpendicular drop method using autopeaks.m function
    ap=autopeaks(xx,yy,1);
    AreaAP(trial)=ap(5);
    
    % Curve fitting with exponentially-broadened Gaussian (shown in
    % lower-left panel)
    NumPeaks=1;
    PeakShape=31; % Exponential broadened Gaussian (variable time constant)
    Pcfeg=peakfit([xx;yy],0,0,NumPeaks,PeakShape,0,0,0,0,0,0);
    ge=Pcfeg(3).*gaussian(xx,Pcfeg(2),Pcfeg(4)); % Original Gaussian recovered
    % figure(1);subplot(2,2,3);plot(xx,yy,'.',xx,ge,'r')
    AreaCfeg(trial)=Pcfeg(5)./(xx(2)-xx(1));
    
    % Curve fitting with two overlapping Gaussians (shown in lower-right panel)
    NumPeaks=2;
    PeakShape=1; % simple symmetrical Gaussian
    Pcd2g=peakfit([xx;yy],0,0,NumPeaks,PeakShape,0,0,0,0,0,0);
    g2a=Pcd2g(1,3).*gaussian(xx,Pcd2g(1,2),Pcd2g(1,4)); % First gaussian
    g2b=Pcd2g(2,3).*gaussian(xx,Pcd2g(2,2),Pcd2g(2,4)); % Second gaussian
    g2=g2a+g2b; % Sum of the two Gaussians
    Areacd2g(trial)=sum(g2)./(xx(2)-xx(1));
   
end
AreaTm=rmnan(AreaTm); % Remove NaNs from this vector
disp('Asymmetrical Area Measurement Precision Test     Standard deviation  ')
disp(['Total area (sum of y''s divided by delta x).......... ' num2str(std(TotalArea))  ])
disp(['Gaussian estimation method (findpeaks.m)............ ' num2str(std(AreaGge)) ] )
disp(['Triangulation method (findpeaksT.m)................. ' num2str(std(AreaTm)) ] )
disp(['Perpendicular drop method (autopeaks.m)............. ' num2str(std(AreaAP))  ] )
disp(['Exponentially-broadened Gaussian ................... ' num2str(std(AreaCfeg)) ] ) 
disp(['Sum of two overlapping symmetrical Gaussians........ ' num2str(std(Areacd2g))  ] )

plot(xx,yy,'.')