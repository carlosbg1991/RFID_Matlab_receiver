function G=NumAT(m,threshold)
% "Numbers Above Threshold": counts the number of adjacent elements in the
% vector 'm' that are greater than or equal to 'threshold'. Returns a
% matrix listing each group of adjacent values, their starting index, the
% number of elements in that group, and the sum of that group, and the
% mean. Example:
% m=[6 5 7 1 2 1 2 1 5 6 5 5 3 2 3 2 1 0 0 1 0 1 2 1 2 1 5 5 7 3 2 3 2 1 1 5 5];
% threshold=4;
% G=NumAT(m,threshold);
% disp('    group     index     count     sum        mean')
% disp(G)
% Result:
%    group     index     count     sum        mean
%     1.0000    1.0000    3.0000   18.0000    6.0000
%     2.0000    9.0000    4.0000   21.0000    5.2500
%     3.0000   27.0000    3.0000   17.0000    5.6667
%     4.0000   36.0000    2.0000   10.0000    5.0000
%
% Tom O'Haver, December 2016
mm=[m 0];
groupcounter=1;
count=0;
sum=0;
clear G
for index=1:length(mm),
    if mm(index)>threshold,
        abovethreshold=1;
        count=count+1;
        sum=sum+mm(index); 
    else
        if abovethreshold,
            G(groupcounter,:)=[groupcounter index-count count sum sum./count];
            groupcounter=groupcounter+1;
            abovethreshold=0;
        end
        count=0;
        sum=0;
    end
end