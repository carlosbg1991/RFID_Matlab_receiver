function y=sinc(x)
% Sinc function, as described on 
% https://www.mathworks.com/help/signal/gs/the-sinc-function.html
y=sin(pi.*x)./(pi.*x);