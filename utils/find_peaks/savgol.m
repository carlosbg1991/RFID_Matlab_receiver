function [y_hat] = savgol(y,SGpar)
%SAVGOL Savitsky-Golay smoothing and differentiation
% Inputs are the matrix of row vectors to be smoothed (y),
% and the optional variables specifying the number of points in
% filter (width), the order of the polynomial (order), and the
% derivative (deriv). The output is the matrix of smoothed
% and differentiated row vectors (y_hat). If number of points,
% polynomial order and derivative are not specified,
% they are set to 15, 2 and 0, respectively. I/O format is:
% y_hat = savgol(y,width,order,deriv);
%
% Example: if y is a 5 by 100 matrix
%          then savgol(y,11,3,1) gives the 5 by 100 matrix of
%          first-derivative row vectors resulting from a 11-
%          point cubic Savitzky-Golay smooth of each row of y

% Sijmen de Jong
% Unilever Research Laboratorium Vlaardingen
% Feb 1993
% Modified by Barry M. Wise
% May 1994, May 1996
%
% Fast version of savgol routine by
% Hans Boelens Feb 1999
%   ==========================================================================
%   Copyright 2005 Biosystems Data Analysis Group ; Universiteit van Amsterdam
%   ==========================================================================
y = y';
width   = SGpar(1);
order   = SGpar(2);
deriv   = SGpar(3);

[m,n]   = size(y);
y_hat   = zeros(m,n);

% In case of input error(s) set to reasonable values
w = max( 3, 1+2*round((width-1)/2) );
if w ~= width
  s = sprintf('Width changed to %g',w);
  disp('  '), disp('Width musth be >= 3 and odd'), disp(s)
end
o = min([max(0,round(order)),5,w-1]);
if o ~= order
  s = sprintf('Order changed to %g',o); disp('  ')
  disp('Order must be <= width -1 and <= 5'), disp(s)
end
d = min(max(0,round(deriv)),o);
if d ~= deriv
  s = sprintf('Derivative changed to %g',d); disp('  ')
  disp('Deriviative must be <= order'), disp(s)
end
p = (w-1)/2;

% Calculate design matrix and pseudo inverse
x = ((-p:p)'*ones(1,1+o)).^(ones(size(1:w))'*(0:o));
weights = ((x'*x)\x')';

% Smoothing and derivative for bulk of the data
ht                  = prod(1:d) * weights(:,d+1);    % Impuls response of SG-filter.
lht                 = length(ht);
temp                = filter(-ht,1.,y');
if ( rem(deriv,2) == 0 )
    y_hat(:,p+2:n-p-1)  = -temp(lht+1:n-2*p+lht-2,:)';
else
    y_hat(:,p+2:n-p-1)  = temp(lht+1:n-2*p+lht-2,:)';
end

% Smoothing and derivative for tails
weights = [y(:,1:w); y(:,n-w+1:n)]*weights;  % full polynomial model
for i=1:d
  weights = weights(:,2:o+2-i)*diag(1:o+1-i); % or its d'th derivative
end
y_hat(:,1:p+1) = weights(1:m,:)*x(1:p+1,1:1+o-d)'; % fitting the tails
y_hat(:,n-p:n) = weights(m+(1:m),:)*x(p+1:w,1:1+o-d)';

y_hat = y_hat';