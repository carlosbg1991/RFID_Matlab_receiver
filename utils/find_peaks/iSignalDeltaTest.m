% Demonstration of the power spectrum of smoothing and differentation
% functions of iSignal, obtained by applying them to a delta function.
disp('Click on the Figure, then use the A and Z keys to increase or ')
disp('decrease the smooth width, the S key to change the smooth type,')
disp('and the D key to change the derivative order.')
x=1:10000;
delta=zeros(size(x));
delta(round(size(delta)./2))=1;delta(1)=0;
isignal([x;delta],5000,10000,1,1000,0,0,0,0,0,0,0,1);