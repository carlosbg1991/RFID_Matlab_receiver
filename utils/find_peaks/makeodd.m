function ok = makeodd(k)
% Returns the nearest odd integer of each elemenr of vector k
% Example:
% makeodd([1.1 2 3 4.8 5 6 7.7 8 9])
% ans =
%      1     3     3     5     5     7     9     9     9
ok=zeros(size(k));
k=round(k);
for n=1:length(k)
    hk=k(n)./2;
    if hk==round(hk)
        ok(n)=k(n)+1;
    else
        ok(n)=k(n);
    end
end