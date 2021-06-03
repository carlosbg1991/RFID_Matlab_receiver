function g = exppulse(x,t1,t2)
% Exponential pulse of the form
% Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
e=(x-t1)./t2;
p = 4*exp(-e).*(1-exp(-e));
p=p .* [p>0];
g = p';