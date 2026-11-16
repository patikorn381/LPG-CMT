function n = neff_core(lam, a1, n1, n2)
% Effective index of the fundamental core mode LP01 (two-layer model).
f  = @(nn) disp2_core(nn, lam, a1, n1, n2);
ns = linspace(n1-1e-9, n2+1e-9, 4000);   % scan from the top down
v  = f(ns);
i  = find(sign(v(1:end-1)) ~= sign(v(2:end)), 1, 'first');
n  = bisect(f, ns(i+1), ns(i));
end
