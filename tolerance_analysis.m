% TOLERANCE_ANALYSIS  One-at-a-time sensitivity of resonant wavelength,
% band depth and bandwidth to fabrication tolerances:
%   dn +-5 %, L +-5 %, Lambda +-1 um, a1 +-0.05 um   (n3 = 1.00)
% A row of NaN means the band left the 780-1220 nm window or vanished
% (for nu = 15: the grating period moved past the turning point).
% Output: results/tolerance.csv with
%   param (1 = dn, 2 = L, 3 = Lambda, 4 = a1), change, nu, lambda_r, d_lambda_r,
%   Tmin_dB, d_Tmin_dB, FWHM_nm
clear; P = lpg_params(); outdir = 'results'; if ~exist(outdir,'dir'), mkdir(outdir); end
t0 = tic;
names = {'dn','L','Lambda','a1'};
cases = [1 -0.05; 1 0.05; 2 -0.05; 2 0.05; 3 -1; 3 1; 4 -0.05; 4 0.05];
[~, ~, r0] = lpg_spectrum(P, P.M, 0.78, 1.22);
out = [];
for c = 1:size(cases,1)
    Q = P; p = cases(c,1); v = cases(c,2);
    switch p
        case 1, Q.dn  = P.dn*(1+v);
        case 2, Q.L   = P.L*(1+v);
        case 3, Q.Lam = P.Lam + v;
        case 4, Q.a1  = P.a1 + v;
    end
    [~, ~, r] = lpg_spectrum(Q, Q.M, 0.78, 1.22);
    for nu = [9 11 13 15]
        i0 = find(r0(:,1) == nu, 1, 'first'); i1 = find(r(:,1) == nu, 1, 'first');
        if isempty(i0), continue; end
        if isempty(i1)   % band moved out of 780-1220 nm or vanished (period beyond turning point)
            out(end+1,:) = [p v nu NaN NaN NaN NaN NaN]; %#ok<AGROW>
            continue;
        end
        out(end+1,:) = [p v nu r(i1,2) r(i1,2)-r0(i0,2) r(i1,5) r(i1,5)-r0(i0,5) r(i1,6)]; %#ok<AGROW>
    end
    fprintf('[tol] case %d/%d (%s %+g) done\n', c, size(cases,1), names{p}, v);
end
fprintf('[tol] param change   nu   lambda_r   d_lambda   Tmin(dB)  d_Tmin   FWHM\n');
for i = 1:size(out,1)
    fprintf('      %-6s %+6.2f  %3d  %8.1f  %+8.2f  %8.2f  %+7.2f  %6.1f\n', names{out(i,1)}, out(i,2), out(i,3:end));
end
write_csv(fullfile(outdir,'tolerance.csv'), 'param,change,nu,lambda_r_nm,d_lambda_r_nm,Tmin_dB,d_Tmin_dB,FWHM_nm', out);
fprintf('tolerance_analysis done in %.1f s\n', toc(t0));
