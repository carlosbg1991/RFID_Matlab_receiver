function RS=RSquared(polycoeff,x,y)
    % Compute the correlation coefficient and R-Squared
    if IsOctave,
        cc=corr(polyval(polycoeff,x'),y');
        RS=cc.^2;
    else
        cc=corrcoef(polyval(polycoeff,x),y);
        RS=cc(2).^2;
    end %   if IsOctave,