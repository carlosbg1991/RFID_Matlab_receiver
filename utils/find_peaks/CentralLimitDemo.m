% Script that demonstrates the central limit theorem by showing the
% probability distribution of random variables composed of varying
% numbers of non-normally-distributed components.
% Because real-world quantities are often the balanced sum of many 
% unobserved random events, the central limit theorem provides a partial 
% explanation for the prevalence of the normal probability distribution.
% In the expressions below, the rand functions provides a uniformly-
% distributed variable. You can include additional functions, such as
% sqrt(rand), sin(rand), rand^2, log(rand), to obtain other distributions.
% No matter what the distribution of the single variable, but the
% time you combine four of them, the resulting distribution is normal.
NumPoints=100000; % The length of the random variable arrays.
NumBins=50; % The numnber of bins in the histograms
subplot(2,2,1)
SingleVariable=rand(NumPoints,1);
hist(SingleVariable,NumBins)
title('Single random variable')
xlabel(['STD = ' num2str(std(SingleVariable)) ])
subplot(2,2,2) 
TwoVariables=rand(NumPoints,1)-rand(NumPoints,1);
hist(TwoVariables,NumBins)
title('Two random variables')
xlabel(['STD = ' num2str(std(TwoVariables))])
subplot(2,2,3) 
FourVariables=rand(NumPoints,1)-rand(NumPoints,1)-rand(NumPoints,1)+rand(NumPoints,1);
hist(FourVariables,NumBins)
StandardDeviation4=std(FourVariables);
xlabel(['STD = ' num2str(StandardDeviation4)])
title('Four random variables')
% The 4th plot shows a normal distrubution obtained by using the
% randn function, for comparison.
subplot(2,2,4) 
RANDNvariable=randn(NumPoints,1);
hist(RANDNvariable,NumBins)
title('Normal distribution')
StandardDeviationRANDN=std(RANDNvariable);
xlabel(['STD = ' num2str(StandardDeviationRANDN)])
