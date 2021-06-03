function sd=stdev(a)
% Octave and Matlab compatible standard deviation function
sa=size(a);
reshape(a,1,length(a));
if sa(1)>sa(2),
  sd=std(a);
else
  sd=(std(a'));
end;