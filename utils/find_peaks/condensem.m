function sm=condensem(M,n)
% Condense Matrix.  Condenses the rows of matrix M by a factor of n, where n 
% is a non-zero positive integer.  Produces an approximate version of matrix M,  
% with each group of n adjacent values in each row of M replaced by its average.
% The number of columns is unchanged. Use for reducing the length and
% processing time of sets of over-sampled signals stored in a matrix. Also
% works for vectors, but condense(y,n) is faster for vectors.
% Example: condensem([1 2 3 4 ; 4 5 6 7],2) yields [1.5 3.5 ; 4.5 6.5]
D=size(M);
n=round(n);
m=floor(length(M)/n);
sm=zeros(D(1),m);
for k=1:D(1),
  y=M(k,:);
  if n > 1
      sm(k,:)=mean(reshape(y(1:n*m),n,m));
  else
      sm(k,:)=y;
  end
end