% Matlab/Octave code example

function derivdemo1

% Create an independent variable ("x-axis")
increment=1;
x=[1:increment:300];

% Create a sigmoidal ("S-shaped") signal as a function of x, centered at x=150
y=sigmoid(x,150,10);

% Plot the signal and its 1st, 2nd, and 3rd derivatives in the four quadrants of the plot window
subplot(2,2,1)
plot(x,y)
subplot(2,2,2)
plot(x,deriv(y))
subplot(2,2,3)
plot(x,deriv2(y))
subplot(2,2,4)
plot(x,deriv3(y))

% Internal sub-functions

function d=deriv(a)
% First derivative of vector using 2-point central difference.
% Example: deriv([1 1 1 2 3 4]) yeilds [0 0 .5 1 1 1]
%  T. C. O'Haver, 1988.
n=length(a);
d=zeros(size(a));
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1;
  d(j)=(a(j+1)-a(j-1)) ./ 2;
end

function d=deriv2(a)
% Second derivative of vector using 3-point central difference.
%  T. C. O'Haver, 2006.
n=length(a);
for j = 2:n-1;
  d(j)=a(j+1) - 2.*a(j) + a(j-1);
end
d(1)=d(2);
d(n)=d(n-1);

function d3=deriv3(a)
% Third derivative of vector a
%  T. C. O'Haver, 2008.
n=length(a);
for j = 3:n-2;
  d3(j)=a(j+2) - 2.*a(j+1) + 2.*a(j-1) - a(j-2);
end
d3(1:2)=d3(3);
d3(n-1:n)=d3(n-2);

function g=sigmoid(x,t1,t2)
g=1./(1 + exp((t1 - x)./t2));