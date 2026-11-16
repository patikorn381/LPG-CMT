function [f, A, B] = disp3_clad(X, lam, a1, a2, n1, n2, n3)
% Three-layer (core/cladding/surrounding) scalar LP0m dispersion function
% for cladding modes (n3 < neff < n2). X = u2*a2 is the transverse
% parameter in the cladding, which spaces the roots roughly by pi.
k  = 2*pi/lam;
u2 = X/a2;
n  = sqrt(n2^2 - (u2/k).^2);
u1 = k*sqrt(n1^2 - n.^2);
w  = k*sqrt(n.^2 - n3^2);
% field: core J0(u1 r); cladding A J0(u2 r) + B Y0(u2 r); outside C K0(w r)
J0c = besselj(0,u1*a1); J1c = besselj(1,u1*a1);
J0a = besselj(0,u2*a1); J1a = besselj(1,u2*a1);
Y0a = bessely(0,u2*a1); Y1a = bessely(1,u2*a1);
d = u2.*(J0a.*Y1a - Y0a.*J1a);
A = (J0c.*u2.*Y1a - Y0a.*u1.*J1c)./d;
B = (J0a.*u1.*J1c - u2.*J1a.*J0c)./d;
psi  = A.*besselj(0,u2*a2) + B.*bessely(0,u2*a2);
dpsi = -u2.*(A.*besselj(1,u2*a2) + B.*bessely(1,u2*a2));
f = dpsi.*besselk(0,w*a2,1) + w.*besselk(1,w*a2,1).*psi;
end
