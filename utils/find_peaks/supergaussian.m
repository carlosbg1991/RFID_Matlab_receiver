function SG=supergaussian(x,position,width,n)
SG=exp(-2.*((x-position)./width).^n);
