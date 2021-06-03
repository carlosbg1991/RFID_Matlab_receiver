 A=1;
s=100;
lambda=.005; % first exponential modification
tau=1/lambda;
mu=500;
deltat=5; % x-axis increment
t=0:deltat:3500;
G=A.*gaussian(t,mu,s*(2*sqrt(2*log(2))));
EMG=A.*expgaussian2(t,mu,s,tau);
FDG=derivxy(t,G);
GMFD=G-tau.*derivxy(t,EMG);
plot(t,G,t,FDG,t,EMG,'g',t,GMFD,'r.')

