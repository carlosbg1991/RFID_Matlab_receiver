function a=rmnan(a)
% Removes NaNs and Infs from vectors, replacing with neighboring real numbers.
% Example:
%  >> v=[1 2 3 4 Inf 6 7 Inf  9];
%  >> rmnan(v)
%  ans =
%     1     2     3     4     4     6     7     7     9
% (c) Tom O'Haver, toh@umd.edu. Version 2, December 2017
la=length(a);
changes=0;
if isnan(a(1)) || isinf(a(1)),a(1)=0;end
for point=1:la
    if isnan(a(point)) || isinf(a(point))
        a(point)=a(point-1);
        changes=changes+1;
    end
end
if changes>0
    disp([num2str(changes) ' changes made.' ])
end