% Demonstration of the running histogram of an expanding sample of randome
% numbers drawn from a population with a mean and a standard deviation of
% 1.000. Two conclusions: (1) distribution approaches a Gaussian only for
% very larger samples, and (2) for small samples, such as 30 or below) you
% get very poor estimates of the true population mean and the standard
% deviation.
maxcount=3; % Maximum sample size (try numbers from 3 to 100)
sumy=0;
cnt=0;
NumberOfBins=20;
clear A
figure(1)
clf
for cnt=1:maxcount % cnt = Sample size (1 to maxcount)
    y = 1+randn(1); % Random number with mean and standard deviation of 1.
    A(cnt)=y;
    histogram(A,NumberOfBins) % Plots the histogram
    xlabel('Value of individual samples')
    ylabel('Number within each bin')
    title(['Count = ' num2str(cnt) '       Average = ' num2str(mean(A))  '       Standard deviation = ' num2str(std(A)) ])
    drawnow
end
% At the end, fit the histogram to a Gaussian model
figure(2)
[N,X]=hist(A,20); % Fit histogram to Gaussian
[F,E]=peakfit([X;N]);
subplot(2,1,1)
title('Fitting a Gaussian model to the histogram data')
R2=E(2); % Coefficent of determination for Gaussian fit.
disp('Results of fitting Gusssian model to the histogram:')
disp(['Peak center = ' num2str(F(2)) ])
CalculatedSigma=F(4)/(2*sqrt(2*log(2)));
disp(['Standard deviation of distribution = ' num2str(CalculatedSigma) ])
disp(['Coefficient of determination = ' num2str(E(2)) ])
%  figure(1)
%  title(['Count = ' num2str(cnt) '       Average = ' num2str(mean(A)) '      Gaussian R2 = ' num2str(R2) ])