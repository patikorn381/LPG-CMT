function w = band_fwhm(lam, Tm, q)
% Full width (nm) of an attenuation band at half of its depth in linear
% transmission, i.e. where Tm = (1 + Tmin)/2. Tm is the single-mode
% transmission and q the grid index at phase matching. Returns NaN when the
% band is cut by the edge of the wavelength window.
level = 0.5*(1 + Tm(q));
iL = q; while iL > 1 && Tm(iL) < level, iL = iL - 1; end
iR = q; while iR < numel(Tm) && Tm(iR) < level, iR = iR + 1; end
if iL == 1 || iR == numel(Tm), w = NaN; return; end
lL = interp1(Tm([iL iL+1]), lam([iL iL+1]), level);
lR = interp1(Tm([iR-1 iR]), lam([iR-1 iR]), level);
w = (lR - lL)*1e3;
end
