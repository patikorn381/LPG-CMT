function n = sellmeier_silica(lam)
% Refractive index of fused silica, Malitson (1965) Sellmeier equation.
% lam in micrometres (valid 0.21-3.71 um).
l2 = lam.^2;
n = sqrt(1 + 0.6961663*l2./(l2 - 0.0684043^2) ...
           + 0.4079426*l2./(l2 - 0.1162414^2) ...
           + 0.8974794*l2./(l2 - 9.896161^2));
end
