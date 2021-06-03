x=0:.01:20;y=zeros(size(x));
y(900:1100)=1;                % Create a rectangular function y, 200 points wide
y=y+.01.*randn(size(y));      % Noise added before the convolution
c=exp(-(1:length(y))./150);    % exponential trailing convolution function, c
yc=conv(y,c,'full')./sum(c);  % Create exponential trailing rectangular function, yc
% yc=yc+.01.*randn(size(yc)); % Noise added after the convolution
ydc=deconv(yc,c).*sum(c);     % Attempt to recover y by deconvoluting c from yc
% Plot all the steps
subplot(2,2,1); 
plot(x,y); 
title('original y');
subplot(2,2,2); 
plot(x,c);
title('c'); 
subplot(2,2,3); 
plot(x,yc(1:2001)); 
title('yc'); 
subplot(2,2,4); 
plot(x,ydc);
title('recovered y')