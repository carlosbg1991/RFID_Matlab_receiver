function y=quadslope(x,boa,coa) % normalized quadratic
 y=(x.^2+(boa).*x+coa);
 y=y./max(y);