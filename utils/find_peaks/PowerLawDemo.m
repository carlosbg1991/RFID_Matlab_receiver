% Simple demonstration of the power transform method. Plots noisy Gaussian
% raised to the power p, normalized to 1.0. Shows that as the power
% increases, peak width decreases and noise is reduced on the baseline but
% increased on the peak maximum.
figure(1)
clf
width=5;
x=0:.005:15;
noise=randn(size(x));

for p=1:2:8
    subplot(2,2,.5+p/2)
    g=100.*gaussian(x,7,width)+noise;
    y=g.^p; % Raise to power p
    y=y./max(y); % Normalized to 1.0
    plot(x,y)
    w(p)=halfwidth(x,y);
    text(0.1,.9,['power=' num2str(p) ])
    text(0.1,.7,['width=' num2str(w(p)) ])
    title(['Normalized y^' num2str(p)'])
end

