function kap = coupling(lam, nco, ncl, a1, a2, n1, n2, n3, dn, nr)
% Coupling coefficient between the LP01 core mode and each cladding mode
% (LP02, LP03, ...), eq. (3):
%
%   kappa = (pi*dn/lam) * Int_0^a1 Eco*Ecl r dr
%                        / sqrt( Int_0^a2 Eco^2 r dr * Int_0^a2 Ecl^2 r dr )
%
% Scalar, weak-guidance form. dn is the AMPLITUDE of the index modulation,
% uniform over the core (the photosensitive region). The DC (self-coupling)
% index change is neglected. The evanescent tail outside a2 is neglected in
% the normalisation of the cladding mode.
%
% The integrals use the trapezoidal rule on two separate grids, [0,a1] and
% [a1,a2], so that the core boundary is a grid node (a single grid over
% [0,a2] truncates the core integral by up to one step and makes kappa
% depend on the grid at the 0.1 % level). nr = number of points of the
% cladding grid (default 12,001); the core grid has nr/6 + 1 points.
if nargin < 10 || isempty(nr), nr = 12001; end
k   = 2*pi/lam;
rc  = linspace(0,  a1, round((nr-1)/6) + 1);     % core grid
rl  = linspace(a1, a2, nr);                      % cladding grid
u   = k*sqrt(n1^2 - nco^2); w = k*sqrt(nco^2 - n2^2);
Ec  = besselj(0,u*rc)/besselj(0,u*a1);
El  = besselk(0,w*rl,1)./besselk(0,w*a1,1).*exp(-w*(rl-a1));
Nco = trapz(rc, Ec.^2.*rc) + trapz(rl, El.^2.*rl);
kap = nan(size(ncl));
for m = 1:numel(ncl)
    if isnan(ncl(m)), continue; end
    X = a2*k*sqrt(n2^2 - ncl(m)^2);
    [~, A, B] = disp3_clad(X, lam, a1, a2, n1, n2, n3);
    u1 = k*sqrt(n1^2 - ncl(m)^2); u2 = X/a2;
    Fc = besselj(0,u1*rc);
    Fl = A*besselj(0,u2*rl) + B*bessely(0,u2*rl);
    Ncl = trapz(rc, Fc.^2.*rc) + trapz(rl, Fl.^2.*rl);
    kap(m) = abs(pi*dn/lam * trapz(rc, Ec.*Fc.*rc) / sqrt(Nco*Ncl));
end
end
