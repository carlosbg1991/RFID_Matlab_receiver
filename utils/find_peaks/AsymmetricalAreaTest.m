% Asymmetrical Area Measurement Test. Test of accuracy of peak area
% measurement methods for an asymmetrical peak,comparing (A) Gaussian
% estimation,(B) triangulation, (C) perpendicular drop method, and curve
% fitting by (D) exponentially braodened Gaussian, and (E) two overlapping
% Gaussians. Must have the following functions in the Matlab/Octave path:
% gaussian.m, expgaussian.m, findpeaksplot.m, findpeaksTplot.m,
% autopeaks.m, and peakfit.m.
clf
disp(' ')
xx=1:100;
% Uncomment one of these lines,13,15,17,or 19, to activate that shape.
% Alternative shape: Two closely overlapping Guassians
yy=2.*gaussian(xx,max(xx)/2.4,max(xx)/4)+gaussian(xx,max(xx)/2,max(xx)/3);
% Alternative shape: exponentially-broadened Gaussian
% yy=2.*expgaussian(xx,max(xx)/2.4,max(xx)/4,-7)';
% Alternative shape: triangle
% yy=2.*triangle(xx,max(xx)/2.4,max(xx)/4);
% Alternative shape: single Gaussian
% yy=2.*gaussian(xx,max(xx)/2.4,max(xx)/4);

noise=0.01; % Increase noise to test for noise robustness
yy=yy+noise.*randn(size(xx));
TotalArea=sum(yy)./(xx(2)-xx(1));

figure(1);subplot(2,2,1)
Pge=findpeaksplot(xx,yy,.01,1,3,max(xx)/4,1);
AreaGge=Pge(5)./(xx(2)-xx(1));
title('Gaussian estimation method (findpeaks.m)')
hold on
subplot(2,2,1)
ge=Pge(3).*gaussian(xx,Pge(2),Pge(4));
plot(xx,ge,'r')
hold off

% Triangulation method (shown in upper-right panel)
subplot(2,2,2)
Ptm=findpeaksTplot(xx,yy,.01,1,3,5,1);
AreaTm=Ptm(5)./(xx(2)-xx(1));
title('Triangulation method (findpeaksT.m)')

% Perpendicular drop method using autopeaks.m function
ap=autopeaks(xx,yy,1);
AreaAP=ap(5);

% figure(2) % Curve fitting with simple Gaussian model
% Pcfg=peakfit([xx;yy]);

% Curve fitting with exponentially-broadened Gaussian (shown in lower-left panel)
figure(3)
NumPeaks=1;
PeakShape=31; % Exponential broadened Gaussian (variable time constant)
Pcfeg=peakfit([xx;yy],0,0,NumPeaks,PeakShape);
ge=Pcfeg(3).*gaussian(xx,Pcfeg(2),Pcfeg(4)); % Original Gaussian recovered
figure(1);subplot(2,2,3);plot(xx,yy,'.',xx,ge,'r')
AreaCfeg=Pcfeg(5)./(xx(2)-xx(1));
title('Symmetrical Gaussian (red) broadened exponentially')

% Curve fitting with two overlapping Gaussians (shown in lower-right panel)
figure(4)
NumPeaks=2;
PeakShape=1; % simple symmetrical Gaussian
Pcd2g=peakfit([xx;yy],0,0,NumPeaks,PeakShape);
g2a=Pcd2g(1,3).*gaussian(xx,Pcd2g(1,2),Pcd2g(1,4)); % First gaussian
g2b=Pcd2g(2,3).*gaussian(xx,Pcd2g(2,2),Pcd2g(2,4)); % Second gaussian
g2=g2a+g2b; % Sum of the two Gaussians
figure(1);subplot(2,2,4);plot(xx,yy,'.',xx,g2a,'g',xx,g2b,'g',xx,g2,'r')
Areacd2g=sum(g2)./(xx(2)-xx(1));
title('Sum of two overlapping symmetrical Gaussians')

disp('Asymmetrical Area Measurement Test          Measured Area   % Error')
disp(['Total area (sum of y''s divided by delta x)    ' num2str(TotalArea) ])
disp(['Gaussian estimation method (findpeaks.m)      ' num2str(AreaGge) '       ' num2str(100.*(TotalArea-AreaGge)./TotalArea) '%'])
disp(['Triangulation method (findpeaksT.m)           ' num2str(AreaTm) '       ' num2str(100.*(TotalArea-AreaTm)./TotalArea) '%'])
disp(['Perpendicular drop method (autopeaks.m)       ' num2str(AreaAP) '       ' num2str(100.*(TotalArea-AreaAP)./TotalArea) '%'])
disp(['Exponentially-broadened Gaussian              ' num2str(AreaCfeg) '       ' num2str(100.*(TotalArea-AreaCfeg)./TotalArea) '%']) 
disp(['Sum of two overlapping symmetrical Gaussians  ' num2str(Areacd2g) '       ' num2str(100.*(TotalArea-Areacd2g)./TotalArea) '%'])