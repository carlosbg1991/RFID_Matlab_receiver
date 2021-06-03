% DerivativeScaling
% Demonstration that the amplitude of the nth derivative of a Gaussian peak
% of height H and width W can be estimated by the empirical equation
% H*(10^(0.027*n^2+n*0.45-0.31)).*W^(-n), where W is the full width at half
% maximum (FWHM) measured in the number of x,y data points.
format compact
x=1:400;
for n=1:6; % n is the derivative order (1 for first derivative, etc)
    for W=10:100, % W is the full width at half maximum of the Gaussian peak
        y=gaussian(x,max(x)./2,W); % Gaussian 
        switch n, % Compute the nth derivative of y
            case 1
                d=deriv1(y); % d is the derivative of y
            case 2
                d=deriv1(deriv1(y));
            case 3
                d=deriv1(deriv1(deriv1(y)));
            case 4
                d=deriv4(y);
             case 5
                d=deriv4(deriv1(y));
             case 6
                d=deriv4(deriv1(deriv1(y)));
        end % switch n
        d(1:n)=0;
        d(length(x)-n:length(x))=0;
        ampd(W)=max(abs(d)); % ampd is the absolute amplitude of d
    end %  for W
    figure(2)
    coeff=plotit((10:100).^-(n),ampd(10:100),1);
    xlabel('W^-^n')
    ylabel('max(abs(d))')
    title('The amplitude of the nth derivative of a peak is inversely proportional to the nth power of its width W')
    slope(n)=coeff(1);
    figure(1)
    plot(x,d)
    xlabel('x')
    ylabel('d')
    pause(1)
end % for n
figure(3)
plotit(1:n,log10(slope),2)
xlabel('Derivative order')
ylabel('Log of the slope of W^(-n) vs max(abs(d))')