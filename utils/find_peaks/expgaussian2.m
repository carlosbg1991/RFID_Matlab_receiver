function emg=expgaussian2(t,mu,s,tau)
if tau~=0,
    lambda=1/tau;
    EMG=s.*lambda.*sqrt(pi/2).*exp(0.5.*(s.*lambda).^2-lambda.*(t-mu)).*erfc((1/sqrt(2)).*(s.*lambda-((t-mu)./s)));
    emg=EMG./max(EMG);
end