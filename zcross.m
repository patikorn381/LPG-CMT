function r = zcross(x, f, x0)
% zero crossing of f(x) nearest to x0 (linear interpolation)
z = find(sign(f(1:end-1)) ~= sign(f(2:end)));
if isempty(z), r = NaN; return; end
xr = x(z) - f(z).*(x(z+1)-x(z))./(f(z+1)-f(z));
[~, i] = min(abs(xr - x0)); r = xr(i);
end
