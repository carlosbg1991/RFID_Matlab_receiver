% Demonstraton of the Law of Large Numbers. Computes the sample standard
% deviation of small samples of numbers drawn from a population of normally
% distributed random numbers with a true population standard deviation of
% 1.00. The sample size (SampleSize) varies from 2 to MaxSampleSize. At each
% SampleSize, NumRepeats samples are measured. The plots shows the standard
% deviation of each sample measured as a dot. The two horizontal red lines
% show an (arbitrary) plus or minus 10% range. Conclusion: you need at
% least a sample size of 10 to get a standard deviation within +-50% and a
% sample size of 100 to get within +-20% of the true standard deviation.
clear
MaxSampleSize=100;
NumRepeats=30;
STD=zeros(1,MaxSampleSize);
SampleSize=3;
        
        for repeats=1:NumRepeats
            y = randn(1,SampleSize);
            STD=std(y);
            histogram(STD)
            hold on
        end
% 
%     axis([1 MaxSampleSize 0 2])
%     grid on
%     xlabel('Sample size')
%     ylabel('Standard Deviation')
%     title('Red lines = -+10%.  Blue lines = -+20%')
%     drawnow
    
