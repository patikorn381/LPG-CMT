function x = bisect(f, a, b, tol, maxit)
% BISECT  Bisection root finder (as described in the paper's methodology).
%   f must change sign on [a,b].
if nargin < 4, tol = 1e-12; end
if nargin < 5, maxit = 200; end
fa = f(a);
for it = 1:maxit
    x = 0.5*(a+b); fx = f(x);
    if fx == 0 || 0.5*(b-a) < tol, return; end
    if sign(fx) == sign(fa), a = x; fa = fx; else, b = x; end
end
end
