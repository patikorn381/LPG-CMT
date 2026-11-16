% DISPERSION_CHECK  Effect of silica material dispersion on the results.
% Baseline: constant n1, n2 (paper). Check: n2(lam) from the Sellmeier
% equation with the core-cladding step n1 - n2 kept constant.
% Reports (a) resonances in air, (b) the SRI sweep, (c) for every SRI value
% whether n3 >= n2(lambda_r), i.e. whether the cladding mode would be leaky.
clear; P = lpg_params(); outdir = 'results'; if ~exist(outdir,'dir'), mkdir(outdir); end
t0 = tic; D = P; D.disp = true;
[~, ~, r0] = lpg_spectrum(P, P.M, 0.78, 1.22);
[~, ~, r1] = lpg_spectrum(D, D.M, 0.78, 1.22);
fprintf('[disp] resonances in air (nu, no-dispersion nm, with-dispersion nm, diff nm, Tmin dB no/with)\n');
tab = [];
for j = 1:size(r0,1)
    q = find(r1(:,1) == r0(j,1), 1, 'first'); if isempty(q), continue; end
    tab(end+1,:) = [r0(j,1) r0(j,2) r1(q,2) r1(q,2)-r0(j,2) r0(j,5) r1(q,5)]; %#ok<AGROW>
end
fprintf('       nu=%2d  %8.2f  %8.2f  %+6.2f   %7.2f  %7.2f\n', tab.');
write_csv(fullfile(outdir,'dispersion_resonances.csv'), 'nu,lr_nodisp_nm,lr_disp_nm,diff_nm,Tmin_nodisp_dB,Tmin_disp_dB', tab);
sri = [P.sri P.sri_ext]; modes = P.sri_modes; lamS = linspace(0.78, 1.30, 53);   % includes the extended point 1.45
l0 = sri_sweep(P, sri, modes, lamS)*1e3;
l1 = sri_sweep(D, sri, modes, lamS)*1e3;
s0 = l0 - repmat(l0(1,:), numel(sri), 1); s1 = l1 - repmat(l1(1,:), numel(sri), 1);
S0 = diff(l0)./repmat(diff(sri(:)),1,4); S1 = diff(l1)./repmat(diff(sri(:)),1,4);
% cladding index at the (baseline) resonance and the leaky flag
n2r = sellmeier_silica(l0/1e3); leaky = repmat(sri(:),1,4) >= n2r;
fprintf('[disp] shift (nm) rel. to air, no dispersion | with dispersion (nu 9 11 13 15); * = SRI >= n2(lambda_r)\n');
for s = 1:numel(sri)
    fl = repmat(' ', 1, 4); fl(leaky(s,:)) = '*';
    fprintf('       SRI %.2f: %s | %s  %s\n', sri(s), sprintf('%8.2f', s0(s,:)), sprintf('%8.2f', s1(s,:)), fl);
end
fprintf('[disp] sensitivity 1.43-1.44 (nm/RIU): no disp %s | disp %s\n', sprintf('%8.1f', S0(end-1,:)), sprintf('%8.1f', S1(end-1,:)));
write_csv(fullfile(outdir,'dispersion_sri.csv'), ...
  'SRI,shift_nodisp_nu9,nu11,nu13,nu15,shift_disp_nu9,nu11,nu13,nu15,n2_at_lr_nu9,nu11,nu13,nu15', ...
  [sri(:) s0 s1 n2r]);
write_csv(fullfile(outdir,'dispersion_sensitivity.csv'), 'SRI_mid,S_nodisp_nu9,nu11,nu13,nu15,S_disp_nu9,nu11,nu13,nu15', ...
  [0.5*(sri(1:end-1)'+sri(2:end)') S0 S1]);
fprintf('dispersion_check done in %.1f s\n', toc(t0));
