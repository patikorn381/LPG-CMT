% ADDED_ANALYSES  Additional analyses reported in the paper.
%  (A) DC self-coupling: a uniform average index change dn_dc in the core,
%      added to the effective indices through the core power fractions
%      (confinement.m). dn_dc = 0 reproduces the model of the paper.
%  (B) Band depth, FWHM and figure of merit FOM = |S|/FWHM versus SRI for
%      the modes nu = 9, 11, 13, 15 (model of the paper, dn_dc = 0).
% Outputs: results/dc_resonances.csv, results/dc_sri.csv,
%          results/sri_band_shape.csv, results/fom.csv
clear; P = lpg_params(); outdir = 'results'; if ~exist(outdir,'dir'), mkdir(outdir); end
t0 = tic;
%% (A1) resonances in air for several dn_dc
dcs = [0 0.5e-4 1e-4 2e-4]; out = [];
for d = dcs
    Q = P; Q.dn_dc = d; Q.nc = 90;
    [~, ~, res] = lpg_spectrum(Q, 8, 0.78, 1.70);
    for r = 1:size(res,1)
        out(end+1,:) = [d res(r,1) res(r,2) res(r,4) res(r,5) res(r,6)]; %#ok<AGROW>
    end
    fprintf('[A1] dn_dc = %.1e : %s\n', d, sprintf('nu%d %.1f nm (%.1f dB)  ', res(:,[1 2 5]).'));
end
write_csv(fullfile(outdir,'dc_resonances.csv'), 'dn_dc,nu,lambda_r_nm,kappaL_rad,Tmin_dB,FWHM_nm', out);
% core-mode confinement and cladding-mode confinement at 1 um (for the text)
[n1, n2] = index_model(P, 1.0); nco = neff_core(1.0, P.a1, n1, n2);
ncl = neff_clad(1.0, P.a1, P.a2, n1, n2, P.n3, 8);
[Gco, Gcl] = confinement(1.0, nco, ncl, P.a1, P.a2, n1, n2, P.n3);
fprintf('[A1] Gamma_core(1 um) = %.4f ; Gamma_clad nu=1..15 = %s\n', Gco, sprintf('%.2e ', Gcl));
%% (A2) SRI sweep with dn_dc = 1e-4 and 2e-4 (first zero = lower branch)
lamS = linspace(0.78, 1.70, 93); modes = P.sri_modes; sri = P.sri; dcsri = [];
for d = [0 1e-4 2e-4]
    Q = P; Q.dn_dc = d;
    lr = sri_sweep(Q, sri, modes, lamS);
    sh = (lr - repmat(lr(1,:), numel(sri), 1))*1e3;
    se = diff(lr*1e3)./repmat(diff(sri(:)), 1, numel(modes));
    dcsri = [dcsri; repmat(d, numel(sri), 1) sri(:) lr*1e3 sh]; %#ok<AGROW>
    fprintf('[A2] dn_dc = %.1e : lr(air) %s nm | shift(1.44) %s nm | S(1.00-1.05) %s | S(1.43-1.44) %s nm/RIU\n', d, ...
        sprintf('%8.2f', lr(1,:)*1e3), sprintf('%8.2f', sh(end,:)), sprintf('%8.2f', se(1,:)), sprintf('%9.1f', se(end,:)));
    smid = 0.5*(sri(1:end-1)' + sri(2:end)');
    for j = 1:numel(modes)
        a = abs(se(:,j)); i1 = find(a >= 10, 1, 'first');
        if isempty(i1), lo = NaN; elseif i1 == 1, lo = sri(1);
        else, lo = interp1(log(a([i1-1 i1])), smid([i1-1 i1]), log(10)); end
        fprintf('[A2]    nu = %d : |S| >= 10 nm/RIU from SRI %.3f\n', 2*modes(j)-1, lo);
    end
end
write_csv(fullfile(outdir,'dc_sri.csv'), 'dn_dc,SRI,lr_nu9,lr_nu11,lr_nu13,lr_nu15,sh_nu9,sh_nu11,sh_nu13,sh_nu15', dcsri);
%% (B) band shape and figure of merit versus SRI (dn_dc = 0)
shape = nan(numel(sri), 1 + 3*numel(modes));
for s = 1:numel(sri)
    Q = P; Q.n3 = sri(s);
    [~, ~, res] = lpg_spectrum(Q, 8, 0.78, 1.22);
    shape(s,1) = sri(s);
    for j = 1:numel(modes)
        r = find(res(:,1) == 2*modes(j)-1, 1, 'first');
        if ~isempty(r), shape(s, 1+3*(j-1)+(1:3)) = res(r, [2 5 6]); end
    end
    fprintf('[B] SRI %.2f : %s\n', sri(s), sprintf('%.2f nm %.1f dB %.2f nm | ', shape(s,2:end)));
end
write_csv(fullfile(outdir,'sri_band_shape.csv'), ...
  'SRI,lr9,Tmin9_dB,FWHM9,lr11,Tmin11_dB,FWHM11,lr13,Tmin13_dB,FWHM13,lr15,Tmin15_dB,FWHM15', shape);
lr  = shape(:, 2:3:end); fw = shape(:, 4:3:end);
S   = diff(lr)./repmat(diff(sri(:)), 1, numel(modes));
fwm = 0.5*(fw(1:end-1,:) + fw(2:end,:));
fom = abs(S)./fwm;                                   % 1/RIU
smid = 0.5*(sri(1:end-1)' + sri(2:end)');
write_csv(fullfile(outdir,'fom.csv'), 'SRI_mid,S9,S11,S13,S15,FOM9,FOM11,FOM13,FOM15', [smid S fom]);
fprintf('[B] FOM (1/RIU) at 1.00-1.05: %s\n', sprintf('%.3f ', fom(1,:)));
fprintf('[B] FOM (1/RIU) at 1.30-1.35: %s\n', sprintf('%.3f ', fom(7,:)));
fprintf('[B] FOM (1/RIU) at 1.43-1.44: %s\n', sprintf('%.3f ', fom(end,:)));
fprintf('added_analyses done in %.0f s\n', toc(t0));
