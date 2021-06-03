% Demonstration that a non-Gaussian peak can be fit with multiple Gaussian
% components, and that the total area of the components approaches the area
% under the non-Gaussian peak as the number of components increases.
% Requires peakfit.m and the non-gaussian shape function (e.g. lorentzian,
% alphafunction, expgaussian) in the path.
x=1:300;
% Select one definition for y
% y=lorentzian(x,150,25)';
% y=alphafunction(x,20,10)';
y=expgaussian(x,40,10,-20);
for n=1:5 % n=number of components
    shape=1; % Gaussian fitting components
    NumTrials=30;
    start=[42 11 48 13 56 19 68 31 77 40]; % length must be 2 x maximum n
    [FitResults,GOF]=peakfit([x' y],0,0,n,shape,0,NumTrials,start(1:2*n));
    disp(FitResults);
    disp(GOF);
    fiterror(n)=GOF(1);
    area(n)=sum(FitResults(:,5));
    AreaError(n)=100*(sum(y)-area(n))./sum(y);
    drawnow
end
clf
plot(1:n,fiterror,1:n,AreaError,1:n,area,1:n,sum(y)*ones(5))
xlabel('Number of Gaussian peaks in model')
text(1,fiterror(1),' Percent fitting error (blue)')
text(1,AreaError(1),' Percent error in area estimation (red)')
text(1,area(1),' Sum of areas of components (yellow)')
text(1,sum(y),' True area (cyan)')
title('Number of Gaussian components needed to fit a non-gaussian peak')

