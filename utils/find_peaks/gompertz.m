function y=gompertz(t,Bo,Kh,L)
% function y=Gompertz(t,Bo,Kh,L)
% A Gompertz curve or Gompertz function, named after Benjamin Gompertz, is
% a sigmoid function. It is a type of mathematical model for a time series,
% where growth is slowest at the start and end of a time period. The
% right-hand or future value asymptote of the function is approached much
% more gradually by the curve than the left-hand or lower valued asymptote,
% in contrast to the simple logistic function in which both asymptotes are
% approached by the curve symmetrically. It is a special case of the
% generalized logistic function.
% 
% Example:
%  x=1:.1:10;y=gompertz(x,10,2,3);plot(x,y)

y=Bo*exp(-exp((Kh*exp(1)/Bo)*(L-t) +1));