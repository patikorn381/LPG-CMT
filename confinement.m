function [Gco, Gcl] = confinement(lam, nco, ncl, a1, a2, n1, n2, n3, nr)
% Fraction of the power of the LP01 core mode (Gco) and of each cladding
% mode (Gcl) that lies inside the core (0 <= r <= a1). Used for the DC
% self-coupling term: a uniform average index change dn_dc in the core
% raises the effective index of a mode by approximately dn_dc*G (first-order
% perturbation, scalar weak-guidance form; Erdogan, JLT 15, 1277 (1997)).
% Same fields and grids as coupling.m.
if nargin < 9 || isempty(nr), nr = 12001; end
k   = 2*pi/lam;
rc  = linspace(0,  a1, round((nr-1)/6) + 1);
rl  = linspace(a1, a2, nr);
u   = k*sqrt(n1^2 - nco^2); w = k*sqrt(nco^2 - n2^2);
Ec  = besselj(0,u*rc)/besselj(0,u*a1);
El  = besselk(0,w*rl,1)./besselk(0,w*a1,1).*exp(-w*(rl-a1));
Pin = trapz(rc, Ec.^2.*rc);
Gco = Pin/(Pin + trapz(rl, El.^2.*rl));
Gcl = nan(size(ncl));
for m = 1:numel(ncl)
    if isnan(ncl(m)), continue; end
    X = a2*k*sqrt(n2^2 - ncl(m)^2);
    [~, A, B] = disp3_clad(X, lam, a1, a2, n1, n2, n3);
    u1 = k*sqrt(n1^2 - ncl(m)^2); u2 = X/a2;
    Fc = besselj(0,u1*rc);
    Fl = A*besselj(0,u2*rl) + B*bessely(0,u2*rl);
    Pc = trapz(rc, Fc.^2.*rc);
    Gcl(m) = Pc/(Pc + trapz(rl, Fl.^2.*rl));
end
end
