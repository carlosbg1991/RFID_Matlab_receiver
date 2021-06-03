function [IdentifiedPeaks]=idpeaktable(P,maxerror,Positions,Names)
% Compares the found peak positions in peak table "P" to a database of
% known peaks, in the form of an cell array of known peak maximum positions
% ("Positions") and matching cell array of names ("Names"). If the position
% of a found peak in the signal is closer to one of the known peaks by less
% than the specified maximun error ("maxerror"), that peak is considered a
% match and its peak position, name, error, and amplitude are entered into
% the output cell array "IdentifiedPeaks".
% 
% For information on cell arrays, see
% http://www.mathworks.com/help/matlab/cell-arrays.html
%
% Use "cell2mat" to access numeric elements of IdentifiedPeaks, e.g.
% cell2mat(IdentifiedPeaks(2,1)) returns the position of the first 
% identified peak, cell2mat(IdentifiedPeaks(2,2)) returns its name, etc.
% Version 1, February 14, 2015, by Tom O'Haver (toh@umd.edu)
%
% Example: 
% >> load P
% >> load DataTableCu
% >> idpeaktable(P,.01,Positions,Names)
%
% This example loads the peak table "P", loads the database "DataTable"
% (containing an array of known peak maximum positions "Positions" and
% matching cell array of names "Names"), lists any peak in Cu that matches
% a known peaks in Positions within an error tolerance of 0.01 nm, and
% returns a table of identified peaks.
%
% Related function: idpeaks.m
%
%  Column labels for Named Peaks table 
IdentifiedPeaks(1,1)={'Position'};
IdentifiedPeaks(1,2)={'Name'};
IdentifiedPeaks(1,3)={'Error'};
IdentifiedPeaks(1,4)={'Amplitude'};
p=2;
for n=1:length(P(:,2)), % Look at each peak detected
   % m=index of the cloest match in Positions
   m=val2ind(Positions,P(n,2));
   % PositionError=difference between detected peak and nearest peak in table
   PositionError=abs(P(n,2)-Positions(m));
   if PositionError<maxerror, % Only identify the peaks if the error is less than MaxError
       % Assemble indentified peaks into cell array "IdentifiedPeaks"
       IdentifiedPeaks(p,:)=([P(n,2) Names(m) P(n,2)-Positions(m) P(n,3)]);
       p=p+1;
   end % if PositionError
end  % for n
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);


