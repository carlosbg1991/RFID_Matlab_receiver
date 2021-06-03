function bmd=bimodal(x,std,a,b)
% Random noise with white power spectrum, bimodal distribution with the
% peaks of the two distributions at "a" and "b", each with standard
% deviation "std", equal in length to vector x. 
% Tom O'Haver, 2012
%
% Example 1:
% hist(bimodal([1:1000],.1,-.2,.2),30);
% Example 2:
% model=gaussian([1:1000],500,200);
% noisymodel=model+bimodal(model,.1,-.2,.2);
% plot(noisymodel,'.')
% hist(model,30)
n=1;
bmd=zeros(size(x));
while n<length(x),
    if rand>.5,
        bmd(n)=a+std.*randn;
    else
        bmd(n)=b+std.*randn;
    end % if rand>.5,
    n=n+1;
end % while n<length(xx)-1,