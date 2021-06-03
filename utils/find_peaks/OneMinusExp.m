function g = OneMinusExp(x,pos,wid)
% OneMinusExp(x,pos,wid) = 1-exp(-wid.*(x-pos));
% Example:  x=0:10;y=1-exp(-(.5.*x));plot(x,y)
g = 1-exp(-wid.*(x-pos));