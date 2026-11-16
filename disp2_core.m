function g = disp2_core(n, lam, a1, n1, n2)
% Two-layer (core/cladding) scalar LP0 eigenvalue equation, eq. (1),
% written in pole-free form with exponentially scaled K Bessel functions:
%   u J1(u) K0(w) - w K1(w) J0(u) = 0,
%   u = a1 k sqrt(n1^2 - neff^2),  w = a1 k sqrt(neff^2 - n2^2)
% (equivalent to u J1(u)/J0(u) = w K1(w)/K0(w); weak-guidance approximation)
k = 2*pi/lam;
u = a1*k*sqrt(n1^2 - n.^2);
w = a1*k*sqrt(n.^2 - n2^2);
g = u.*besselj(1,u).*besselk(0,w,1) - w.*besselk(1,w,1).*besselj(0,u);
end
