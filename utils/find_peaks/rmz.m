function a=rmz(a)
% Removes zeros from vectors, replacing with nearest non-zero numbers.
% Example:
%  >> v=[1 2 3 4 0 6 7 0 9];
%  >> rmz(v)
%  ans =
%   1     2     3     4     4     6     7     7     9
% Use this to remove zeros from vectors that will subsequently be used as
% the demonominator of a division.
La=length(a);
changes=0;
if a(1)==0,a(1)=1e-50;end
for point=2:La
    if a(point)==0
        a(point)=a(point-1);
        changes=changes+1;
    end
end
if changes>0
    disp([num2str(changes) ' changes made.' ])
end