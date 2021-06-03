% Demonstraton of the Law of Large Numbers. Computes the sample mean
% (average) of small samples of numbers drawn from a population of normally
% distributed random numbers with a true population mean of 100 and a
% standard deviation of 10. The sample size (SampleSize) varies from 2 to
% MaxSampleSize. At each SampleSize, NumRepeats samples are measured. The
% plots shows the mean of each sample measured as a dot. The two horizontal
% red lines show an (arbitrary) plus or minus 10% range. Conclusion: you
% need at least a sample size of 10 to get a mean within +-1/2 of a
% standard deviation and a sample size of 100 to get within +-1/5 of a
% standard deviation of the true mean.
clear
MaxSampleSize=100;
NumRepeats=10;
Mean=zeros(1,MaxSampleSize);
for SampleSize=1:MaxSampleSize 
   
    if SampleSize>1
        
        for repeats=1:NumRepeats
            y = 100+10.*randn(1,SampleSize);
            Mean=mean(y);
            plot(SampleSize,Mean,'.k',[2 SampleSize],[98 98],'-k',[2 SampleSize],[101 101],'-r',[2 SampleSize],[99 99],'-r',[2 SampleSize],[102 102],'-b',[2 SampleSize],[98 98],'-b')
            hold on
        end
        
    end
    
    axis([1 MaxSampleSize 90 110])
    grid on
    xlabel('Sample size')
    ylabel('Mean (Red lines = -+1%.  Blue lines = -+2%')
   
    title('Running mean of increasing size samples of random numbers.')
    drawnow
    
end
hold off