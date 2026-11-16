function x = xhalf(X)
% half-width x_h of a band of coupling strength X = kappa*L: the smallest
% x > 0 with  X^2/(X^2+x^2) sin^2(sqrt(X^2+x^2)) = sin^2(X)/2
g = @(x) X^2./(X^2+x.^2).*sin(sqrt(X^2+x.^2)).^2 - 0.5*sin(X)^2;
xs = linspace(1e-6, 6, 6001); v = g(xs);
i = find(v(1:end-1) > 0 & v(2:end) <= 0, 1, 'first');
x = fzero(g, xs([i i+1]));
end
