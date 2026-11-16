% VALIDATION_SPECS  Table 3 (rows 1-2) and Table 4 (row 1, Rudolph-Neumann).
% Consistency of the chosen fibre parameters with the PS750 specification:
% NA, cut-off wavelength and mode-field diameter at 780 nm. The MFD is given
% with two definitions:
%   (a) 2 x radius where the exact LP01 intensity falls to 1/e^2 (paper v1)
%   (b) Petermann II, the definition used in datasheets / ITU-T G.650.1
clear; P = lpg_params(); outdir = 'results'; if ~exist(outdir,'dir'), mkdir(outdir); end
NA  = sqrt(P.n1^2 - P.n2^2);
lc  = 2*pi*P.a1*NA/2.405*1e3;                       % nm
rows = [];
for lam = [0.78 1.0]
    [n1, n2] = index_model(P, lam);
    nco = neff_core(lam, P.a1, n1, n2);
    k = 2*pi/lam; u = k*sqrt(n1^2-nco^2); w = k*sqrt(nco^2-n2^2);
    V = k*P.a1*sqrt(n1^2-n2^2);
    r = linspace(0, 40, 40001); ic = r <= P.a1;
    E = zeros(size(r)); dE = E;
    E(ic)  = besselj(0,u*r(ic))/besselj(0,u*P.a1);
    dE(ic) = -u*besselj(1,u*r(ic))/besselj(0,u*P.a1);
    sc = exp(-w*(r(~ic)-P.a1))./besselk(0,w*P.a1,1);
    E(~ic)  = besselk(0,w*r(~ic),1).*sc;
    dE(~ic) = -w*besselk(1,w*r(~ic),1).*sc;
    I = (E/E(1)).^2;
    MFD_1e2 = 2*interp1(flip(I), flip(r), exp(-2));
    MFD_PII = 2*sqrt(2*trapz(r, E.^2.*r)/trapz(r, dE.^2.*r));
    b  = (1.1428 - 0.996/V)^2;                       % Rudolph-Neumann [21]
    nRN = sqrt(n2^2 + b*(n1^2 - n2^2));
    rows(end+1,:) = [lam V nco nRN abs(nco-nRN) MFD_1e2 MFD_PII]; %#ok<AGROW>
end
fprintf('[spec] NA = %.4f (spec %.2f-%.2f), cut-off = %.1f nm (spec %d-%d nm)\n', NA, P.spec_NA, lc, P.spec_cutoff);
for i = 1:size(rows,1)
    fprintf('[spec] lam = %.2f um: V = %.3f, neff = %.7f, neff(RN) = %.7f, |diff| = %.1e, MFD(1/e^2) = %.2f um, MFD(Petermann II) = %.2f um\n', rows(i,:));
end
fprintf('[spec] MFD at 780 nm vs spec %.1f-%.1f um: 1/e^2 %s, Petermann II %s\n', P.spec_MFD780, ...
        ifelse_str(rows(1,6) >= P.spec_MFD780(1) && rows(1,6) <= P.spec_MFD780(2)), ...
        ifelse_str(rows(1,7) >= P.spec_MFD780(1) && rows(1,7) <= P.spec_MFD780(2)));
write_csv(fullfile(outdir,'table3_specs.csv'), 'lambda_um,V,neff,neff_RN,abs_diff,MFD_1e2_um,MFD_PetermannII_um', rows);
write_csv(fullfile(outdir,'table3_NA_cutoff.csv'), 'NA,cutoff_nm', [NA lc]);
