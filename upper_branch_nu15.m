% UPPER_BRANCH_NU15  Second (upper-branch) resonance of nu = 15 at
% Lambda = 400 um versus SRI (dual resonance near the turning point).
clear; P = lpg_params(); outdir = 'results';
lamS = linspace(1.30, 1.60, 61); sri = [1.00 1.33 1.40 1.44]; lr = nan(size(sri));
for s = 1:numel(sri)
    d = nan(size(lamS));
    for i = 1:numel(lamS)
        [n1, n2] = index_model(P, lamS(i));
        nc = neff_core(lamS(i), P.a1, n1, n2);
        nl = neff_clad(lamS(i), P.a1, P.a2, n1, n2, sri(s), 8);
        d(i) = (nc - nl(8))*P.Lam - lamS(i);
    end
    lf = linspace(lamS(1), lamS(end), 20000); df = interp1(lamS, d, lf, 'spline');
    z = find(sign(df(1:end-1)) ~= sign(df(2:end)), 1, 'last');
    if ~isempty(z), lr(s) = lf(z) - df(z)*(lf(z+1)-lf(z))/(df(z+1)-df(z)); end
end
fprintf('[upper nu=15] SRI %s : lambda_r = %s nm ; shift = %s nm\n', sprintf('%.2f ', sri), sprintf('%.2f ', lr*1e3), sprintf('%.2f ', (lr-lr(1))*1e3));
write_csv(fullfile(outdir,'upper_branch_nu15.csv'), 'SRI,lambda_r_nm', [sri(:) lr(:)*1e3]);
