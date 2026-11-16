function n = neff_clad(lam, a1, a2, n1, n2, n3, M)
% Effective indices of the first M circularly symmetric cladding modes of
% the three-layer model. With LP01 denoting the core mode, the j-th cladding
% mode is LP0,j+1 (LP02, LP03, ...). In Erdogan's exact vector model it
% corresponds to the odd-numbered cladding mode nu = 2j - 1, the ones that
% couple strongly to the core mode. Returns NaN when n3 >= n2: the cladding modes are then leaky
% and lie outside the validity of this model.
n = nan(1, M);
if n3 >= n2, return; end
k = 2*pi/lam;
Xmax = k*a2*sqrt(n2^2 - n3^2) - 1e-6;     % neff must stay above n3
f  = @(X) disp3_clad(X, lam, a1, a2, n1, n2, n3);
Xs = linspace(1e-3, min(Xmax, pi*(M+3)), 40*(M+3));
v  = f(Xs);
idx = find(sign(v(1:end-1)) ~= sign(v(2:end)));
for j = 1:min(M, numel(idx))
    X = bisect(f, Xs(idx(j)), Xs(idx(j)+1));
    n(j) = sqrt(n2^2 - (X/(k*a2))^2);
end
end
