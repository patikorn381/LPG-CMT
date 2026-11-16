% COMPARE_SINGH  Table 3 (rows 3-5): comparison with Singh [17].
% Parameter set of [17]: a1 = 4.65 um, a2 = 62.5 um, n1 = 1.465, n2 = 1.46,
% n3 = 1, Lambda = 455 um. dn and L are taken from Table 1 (they do not
% change the resonant wavelengths).
%
% Two comparisons are printed:
%  (a) resonant wavelengths: ALL resonances in 1200-1750 nm are listed.
%      [17] reports three dips at 1378.0, 1439.6 and 1613.1 nm and states
%      that they belong to cladding modes m = 1, 3, 5. The paper (v1)
%      matched them to nu = 3, 5, 7 instead. Both assignments are printed.
%  (b) effective indices at 1310 nm against Table 2 of [17].
% NOTE: the reference numbers below were transcribed from [17]; check them
% against the original article before they are used in the paper.
clear; P = lpg_params(); outdir = 'results'; if ~exist(outdir,'dir'), mkdir(outdir); end
S = P; S.a1 = 4.65; S.n1 = 1.465; S.n2 = 1.46; S.Lam = 455; S.n3 = 1.0; S.disp = false;
[~, ~, res] = lpg_spectrum(S, 12, 1.20, 1.75);
fprintf('[Singh] all resonances in 1200-1750 nm (nu, lambda_r [nm], Tmin [dB]):\n');
fprintf('        nu=%2d  %8.1f nm  %7.2f dB\n', res(:,[1 2 5]).');
ref = [1378.0 1439.6 1613.1];
out = nan(6,5); row = 0;
for assign = 1:2
    if assign == 1, nu_ref = [3 5 7]; else, nu_ref = [1 3 5]; end
    for i = 1:3
        row = row + 1; j = find(res(:,1) == nu_ref(i), 1, 'first');
        if ~isempty(j), out(row,:) = [assign ref(i) nu_ref(i) res(j,2) 100*abs(res(j,2)-ref(i))/ref(i)]; end
    end
end
fprintf('[Singh] assignment A (paper v1: nu = 3, 5, 7)\n');
fprintf('        ref %.1f nm vs nu=%d: %.1f nm, error %.2f %%\n', out(1:3,2:5).');
fprintf('[Singh] assignment B (as stated in [17]: m = 1, 3, 5)\n');
fprintf('        ref %.1f nm vs nu=%d: %.1f nm, error %.2f %%\n', out(4:6,2:5).');
write_csv(fullfile(outdir,'table3_singh.csv'), 'assignment,ref_nm,nu,sim_nm,err_percent', out);
write_csv(fullfile(outdir,'singh_all_resonances.csv'), 'nu,lambda_r_nm,kappa_per_um,kappaL_rad,Tmin_dB,FWHM_nm,FWHM_est_nm', res);
% (b) effective indices at 1310 nm
lam = 1.31; ncl = neff_clad(lam, S.a1, S.a2, S.n1, S.n2, S.n3, 5);
ref3 = [1.459903 1.459743 1.459515 1.458846 1.4588468];   % "three-layer" column of Table 2 in [17]
ref2 = [1.459903 1.459716 1.459155 1.458295 1.457134];    % "two-layer" column of Table 2 in [17]
tab = [(1:2:9).' ncl(:) ref3(:) ncl(:)-ref3(:) ref2(:)];
fprintf('[Singh] neff at 1310 nm: nu | this work (scalar 3-layer) | [17] three-layer | difference | [17] two-layer\n');
fprintf('        %2d | %.6f | %.6f | %+.6f | %.6f\n', tab.');
fprintf('[Singh] core mode at 1310 nm: neff = %.6f\n', neff_core(lam, S.a1, S.n1, S.n2));
write_csv(fullfile(outdir,'singh_neff_1310.csv'), 'nu,neff_this_work,neff_ref_three_layer,diff,neff_ref_two_layer', tab);
