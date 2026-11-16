% CONVERGENCE_CHECK  Numerical convergence of the reported results.
% Every discretisation parameter is doubled and the results are compared
% with the default settings:
%   coarse wavelength grid   45 -> 89 points   (spline nodes for neff, kappa)
%   fine wavelength grid     40,001 -> 80,001  (spectrum, band minima)
%   radial grid              12,001 -> 24,001  (overlap integrals)
%   SRI wavelength grid      53 -> 105 points  (phase-matching residual)
% Output: results/convergence.csv (per band) and convergence_sri.csv.
clear; P = lpg_params(); outdir = 'results'; if ~exist(outdir,'dir'), mkdir(outdir); end
t0 = tic;
R = P; R.nc = 89; R.nf = 80001; R.nr = 24001;
[~, ~, r0] = lpg_spectrum(P, P.M, 0.78, 1.22);
[~, ~, r1] = lpg_spectrum(R, R.M, 0.78, 1.22);
out = [];
for j = 1:size(r0,1)
    q = find(r1(:,1) == r0(j,1), 1, 'first'); if isempty(q), continue; end
    out(end+1,:) = [r0(j,1) r0(j,2) r1(q,2) r1(q,2)-r0(j,2) r0(j,4) r1(q,4) r0(j,5) r1(q,5) r1(q,5)-r0(j,5) r0(j,6) r1(q,6)]; %#ok<AGROW>
end
fprintf('[conv] nu | lambda_r default / refined / diff (nm) | kappa*L default / refined | Tmin default / refined / diff (dB) | FWHM default / refined (nm)\n');
fprintf('       %2d | %9.3f %9.3f %+8.4f | %7.4f %7.4f | %8.2f %8.2f %+6.2f | %6.2f %6.2f\n', out.');
fprintf('[conv] max |d lambda_r| = %.4f nm, max |d kappa*L| = %.2e rad, max |d Tmin| = %.2f dB (nu = %d)\n', ...
        max(abs(out(:,4))), max(abs(out(:,6)-out(:,5))), max(abs(out(:,9))), out(abs(out(:,9)) == max(abs(out(:,9))), 1));
write_csv(fullfile(outdir,'convergence.csv'), ...
  'nu,lr_default_nm,lr_refined_nm,d_lr_nm,kL_default,kL_refined,Tmin_default_dB,Tmin_refined_dB,d_Tmin_dB,FWHM_default_nm,FWHM_refined_nm', out);
sri = P.sri; modes = P.sri_modes;
l0 = sri_sweep(P, sri, modes, linspace(0.78, 1.30, 53))*1e3;
l1 = sri_sweep(P, sri, modes, linspace(0.78, 1.30, 105))*1e3;
fprintf('[conv] SRI sweep: max |d lambda_r| between 53 and 105 grid points = %.4f nm (nu 9, 11, 13, 15: %s)\n', ...
        max(abs(l1(:)-l0(:))), sprintf('%.4f ', max(abs(l1-l0), [], 1)));
write_csv(fullfile(outdir,'convergence_sri.csv'), 'SRI,d_lr_nu9_nm,d_lr_nu11_nm,d_lr_nu13_nm,d_lr_nu15_nm', [sri(:) l1-l0]);
fprintf('convergence_check done in %.1f s\n', toc(t0));
