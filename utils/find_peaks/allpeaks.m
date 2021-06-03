function allpeaks(x,y)
% Looks for data points in x,y that have a y value greater than the adjacent points,
% then prints out the peak number, x value, and y value of those points.
peak=0;
for k = 2:length(x)-1
   if y(k) > y(k-1)
       if y(k) > y(k+1)
           peak = peak + 1;
           disp([round(peak) x(k) y(k)]);
       end
   end  
end