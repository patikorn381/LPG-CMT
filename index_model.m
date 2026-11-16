function [n1, n2] = index_model(P, lam)
% Core and cladding indices at wavelength lam.
%   P.disp = false : constant indices (model used in the original paper)
%   P.disp = true  : cladding follows the silica Sellmeier equation and the
%                    core-cladding index step n1 - n2 is kept constant
if isfield(P, 'disp') && P.disp
    n2 = sellmeier_silica(lam);
    n1 = n2 + (P.n1 - P.n2);
else
    n1 = P.n1; n2 = P.n2;
end
end
