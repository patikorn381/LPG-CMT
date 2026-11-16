function lr = sri_sweep(P, sri, modes, lamS)
% Resonant wavelength (um) of selected cladding modes versus surrounding index.
% modes = cladding-mode indices j (j-th cladding mode = LP0,j+1, nu = 2j - 1).
% For every SRI value the phase-matching residual
%   d(lam) = (neff_co - neff_cl)*Lambda - lam
% is computed on the grid lamS, interpolated with a cubic spline onto 20,000
% points and its FIRST zero crossing is taken (lower branch for nu = 15).
% Points where the cladding mode does not exist (n3 >= n2(lam), leaky
% regime) are excluded; if the resonance falls there the result is NaN.
dndc = pget(P, 'dn_dc', 0);
lr = nan(numel(sri), numel(modes));
for s = 1:numel(sri)
    d = nan(numel(lamS), numel(modes));
    for i = 1:numel(lamS)
        [n1, n2] = index_model(P, lamS(i));
        nc = neff_core(lamS(i), P.a1, n1, n2);
        nl = neff_clad(lamS(i), P.a1, P.a2, n1, n2, sri(s), max(modes));
        if dndc ~= 0   % optional DC self-coupling (P.dn_dc)
            [Gco, Gcl] = confinement(lamS(i), nc, nl, P.a1, P.a2, n1, n2, sri(s));
            nc = nc + dndc*Gco; nl = nl + dndc*Gcl;
        end
        d(i,:) = (nc - nl(modes))*P.Lam - lamS(i);
    end
    for j = 1:numel(modes)
        ok = ~isnan(d(:,j));
        last = find(~ok, 1, 'first') - 1;           % contiguous valid range
        if isempty(last), last = numel(lamS); end
        if last < 4, continue; end
        lf = linspace(lamS(1), lamS(last), 20000);
        df = interp1(lamS(1:last), d(1:last,j), lf, 'spline');
        z = find(sign(df(1:end-1)) ~= sign(df(2:end)), 1, 'first');
        if ~isempty(z)
            lr(s,j) = lf(z) - df(z)*(lf(z+1)-lf(z))/(df(z+1)-df(z));
        end
    end
end
end
