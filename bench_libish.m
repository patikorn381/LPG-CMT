% BENCH_LIBISH  Quantitative external benchmark against the open-access
% experiment of Libish et al., J. Optoelectron. Adv. Mater. 13(5), 491-496 (2011):
% LPGs in hydrogen-loaded SMF-28e (core 1.46145, cladding 1.456, core
% diameter 8.2 um, cladding diameter 125 um), periods 550 and 415 um.
% Measured in air: LPG-1 (550 um) LP02 1451, LP03 1497, LP04 1588 nm;
%                  LPG-2 (415 um) LP02 1254, LP03 1284, LP04 1333, LP05 1423, LP06 1610 nm.
% Measured total blue shift, SRI 1 -> ~1.45: LP04 of LPG-1 21.20 nm, LP06 of LPG-2 105.30 nm.
% Only the parameters stated in that paper are used; dn and L do not change
% the resonant wavelengths. Output: results/bench_libish.csv
clear; P = lpg_params(); outdir = 'results';
S = P; S.a1 = 4.1; S.a2 = 62.5; S.n1 = 1.46145; S.n2 = 1.456; S.n3 = 1.0; S.disp = false;
meas = {550, [1451 1497 1588]; 415, [1254 1284 1333 1423 1610]};
out = [];
for g = 1:2
  S.Lam = meas{g,1}; ref = meas{g,2};
  [~, ~, res] = lpg_spectrum(S, 8, 1.10, 1.80);
  for m = 1:numel(ref)
    j = find(res(:,1) == 2*m-1, 1, 'first');      % lower branch of LP0,m+1
    out(end+1,:) = [S.Lam, m+1, ref(m), res(j,2), res(j,2)-ref(m)]; %#ok<AGROW>
  end
end
% total SRI shifts (lower branch), SRI 1.00 -> 1.45 and 1.454
cases = [550 3; 415 5];  meas_shift = [21.20 105.30];
sh = zeros(2,4);
for c = 1:2
  S.Lam = cases(c,1); lamS = linspace(1.10, 1.80, 141);
  lr = sri_sweep(S, [1.00 1.45 1.454], cases(c,2), lamS)*1e3;
  sh(c,:) = [cases(c,:)+[0 1], lr(2)-lr(1), lr(3)-lr(1)];
end
write_csv(fullfile(outdir,'bench_libish.csv'), 'Lambda_um,LP0m,lambda_meas_nm,lambda_sim_nm,diff_nm', out);
write_csv(fullfile(outdir,'bench_libish_shift.csv'), 'Lambda_um,LP0m,shift_1p45_nm,shift_1p454_nm', sh);
disp(out); disp(sh); disp(meas_shift);
