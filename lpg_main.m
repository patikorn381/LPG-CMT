% LPG_MAIN  Power transmission spectrum of a long-period fibre grating (LPG)
% with coupled-mode theory, following the paper's method:
%   1) LP01 core mode, two-layer eigenvalue equation (bisection)    -> Fig. 1
%   2) LP0m (m >= 2) cladding modes, three-layer model (bisection)  -> Fig. 2
%   3) coupling coefficients from the overlap integral, eq. (3)     -> Fig. 3
%   4) phase-matching curves (resonant wavelength vs. period)       -> Fig. 4
%   5) transmission spectrum 800-1200 nm, eqs. (2) and (4)          -> Fig. 5, Table 2
%   6) surrounding-refractive-index (SRI) sweep                     -> Figs. 6-7, Table 5
% Model: scalar (LP) approximation, uniform grating, independent two-mode
% coupling per cladding mode. Material dispersion is OFF by default
% (P.disp = false) to reproduce the paper; see dispersion_check.m.
%
% Figure conventions follow the LPG literature:
%   Fig. 3  kappa/dn (1/um) vs cladding-mode number
%   Fig. 4  resonant wavelength vs grating period, turning points marked
%   Figs. 5-6  transmission in dB vs wavelength in nm
%   Fig. 7  resonance wavelength shift (nm) vs SRI, linear axes
% The m-th cladding mode is labelled LP0,m+1 (LP01 = core mode) together with
% Erdogan's cladding-mode number nu = 2m - 1.
% Runs in MATLAB and GNU Octave. Units: micrometres.
clear; close all;
DO_PLOT = isempty(getenv('LPG_NO_PLOT'));   % set env LPG_NO_PLOT=1 to skip figures
if exist('OCTAVE_VERSION','builtin') && isempty(available_graphics_toolkits()), DO_PLOT = false; end
% Journal figures carry their description in the caption, not in a title:
% titles are off by default (set SHOW_TITLES = true to draw them).
% Figures are written at 300 dpi.
SHOW_TITLES = false;
if SHOW_TITLES, tt = @(s) title(s); else, tt = @(s) []; end
P = lpg_params(); M = P.M;
outdir = 'results'; if ~exist(outdir,'dir'), mkdir(outdir); end
t_start = tic;
write_csv(fullfile(outdir,'params.csv'), 'a1_um,a2_um,n1,n2,n3,dn,L_um,Lambda_um,M,lam0_um,disp', ...
          [P.a1 P.a2 P.n1 P.n2 P.n3 P.dn P.L P.Lam P.M P.lam0 P.disp]);
% 8 clearly distinct colours (MATLAB's default order repeats after 7)
C8 = [0 0.447 0.741; 0.850 0.325 0.098; 0.929 0.694 0.125; 0.494 0.184 0.556;
      0.466 0.674 0.188; 0.301 0.745 0.933; 0.635 0.078 0.184; 0 0 0];
% Mode names: LP01 is the core mode, so the m-th cladding mode is LP0,m+1
% (LP02, LP03, ...) as in the LPG literature; nu = 2m - 1 is Erdogan's number.
mlabel = @(m) sprintf('LP_{0%s} (\\nu = %d)', lp_index(m+1), 2*m-1);

%% ---- 1) Core mode at the reference wavelength (Fig. 1) ------------------
lam0 = P.lam0; [n1, n2] = index_model(P, lam0);
nco0 = neff_core(lam0, P.a1, n1, n2);
k0 = 2*pi/lam0; u = k0*sqrt(n1^2-nco0^2); w = k0*sqrt(nco0^2-n2^2);
r = linspace(0, 15, 3001); E = zeros(size(r)); ic = r<=P.a1;
E(ic) = besselj(0,u*r(ic))/besselj(0,u*P.a1);
E(~ic) = besselk(0,w*r(~ic),1)./besselk(0,w*P.a1,1).*exp(-w*(r(~ic)-P.a1));
I = (E/max(E)).^2; MFR = interp1(flip(I), flip(r), exp(-2));
fprintf('[1] LP01 at %.2f um: neff = %.7f, mode-field radius (1/e^2) = %.3f um (a1 = %.2f um)\n', lam0, nco0, MFR, P.a1);
write_csv(fullfile(outdir,'core_mode.csv'), 'lambda_um,neff,MFR_1e2_um', [lam0 nco0 MFR]);
if DO_PLOT
figure('visible','off'); plot(r, I, 'LineWidth',1.5); hold on;
plot([P.a1 P.a1],[0 1],'k--'); plot([MFR MFR],[0 1],'r:','LineWidth',1.5);
xlim([0 10]); xlabel('Radial position r (\mum)'); ylabel('Normalized intensity');
tt(sprintf('Fundamental core mode LP_{01}, \\lambda = %.2f \\mum', lam0));
legend('|E|^2','core boundary','1/e^2 radius'); grid on;
print(fullfile(outdir,'fig1_core_mode.png'),'-dpng','-r300');
end

%% ---- 2) Cladding modes + 3) coupling coefficients (Figs. 2-3) ----------
ncl0 = neff_clad(lam0, P.a1, P.a2, n1, n2, P.n3, M);
kap0 = coupling(lam0, nco0, ncl0, P.a1, P.a2, n1, n2, P.n3, P.dn);
fprintf('[2-3] cladding modes at %.2f um (nu, neff, kappa [1/um], kappa/dn [1/um])\n', lam0);
for m = 1:8
    fprintf('      LP0%-2d (nu=%2d): %.8f  %.3e  %.4f\n', m+1, 2*m-1, ncl0(m), kap0(m), kap0(m)/P.dn);
end
write_csv(fullfile(outdir,'coupling_at_lam0.csv'), 'nu,neff_clad,kappa_per_um,kappa_over_dn_per_um', ...
          [2*(1:M)'-1 ncl0(:) kap0(:) kap0(:)/P.dn]);
if DO_PLOT
figure('visible','off'); hold on; rr = linspace(0, P.a2, 4001);
for m = 1:5
    X = P.a2*k0*sqrt(n2^2 - ncl0(m)^2);
    [~, A, B] = disp3_clad(X, lam0, P.a1, P.a2, n1, n2, P.n3);
    u1 = k0*sqrt(n1^2 - ncl0(m)^2); u2 = X/P.a2; Ef = zeros(size(rr)); icr = rr <= P.a1;
    Ef(icr) = besselj(0,u1*rr(icr)); Ef(~icr) = A*besselj(0,u2*rr(~icr)) + B*bessely(0,u2*rr(~icr));
    Ef = Ef/sqrt(trapz(rr, Ef.^2.*rr)); plot(rr, Ef.^2, 'LineWidth', 1.2, 'Color', C8(m,:));
end
yl = ylim; plot([P.a1 P.a1], yl, 'k--'); xlim([0 P.a2]);
xlabel('Radial position r (\mum)'); ylabel('Normalized intensity |E|^2'); grid on; box on;
legend([arrayfun(mlabel, 1:5, 'UniformOutput', false), {'core boundary'}]);
tt(sprintf('Cladding modes at \\lambda = %.2f \\mum', lam0));
print(fullfile(outdir,'fig2_cladding_modes.png'),'-dpng','-r300');
figure('visible','off'); bar(2*(1:M)-1, kap0/P.dn);
xlabel('Cladding mode number \nu (LP_{0m}: \nu = 2m - 3)'); ylabel('\kappa/\Deltan (\mum^{-1})');
tt(sprintf('Normalized coupling coefficients, \\lambda = %.2f \\mum', lam0)); grid on;
print(fullfile(outdir,'fig3_coupling.png'),'-dpng','-r300');
end

%% ---- 4) Phase-matching curves (Fig. 4) ---------------------------------
lamPM = linspace(0.6, 1.7, 56); LamPM = nan(numel(lamPM), 8);
for i = 1:numel(lamPM)
    [n1i, n2i] = index_model(P, lamPM(i));
    nc = neff_core(lamPM(i), P.a1, n1i, n2i);
    nl = neff_clad(lamPM(i), P.a1, P.a2, n1i, n2i, P.n3, 8);
    LamPM(i,:) = lamPM(i)./(nc - nl);
end
% turning points: interior maximum of Lambda(lambda), refined on a spline
tap = nan(8,2);
for m = 1:8
    lf = linspace(lamPM(1), lamPM(end), 5501); Lf = interp1(lamPM, LamPM(:,m), lf, 'spline');
    [Lmax, iT] = max(Lf);
    if iT > 1 && iT < numel(lf) && Lmax <= 800, tap(m,:) = [Lmax lf(iT)*1e3]; end   % inside the plotted range
end
for m = find(~isnan(tap(:,1))).'
    fprintf('[4] turning point of LP0%d (nu = %d): Lambda = %.1f um at lambda = %.0f nm\n', m+1, 2*m-1, tap(m,1), tap(m,2));
end
write_csv(fullfile(outdir,'phase_matching.csv'), 'lambda_nm,Lambda_nu1,nu3,nu5,nu7,nu9,nu11,nu13,nu15', [lamPM(:)*1e3 LamPM]);
write_csv(fullfile(outdir,'turning_points.csv'), 'nu,Lambda_TAP_um,lambda_TAP_nm', [2*(1:8)'-1 tap]);
% both resonances of the highest mode at the working period (dual-resonance)
lf = linspace(lamPM(1), lamPM(end), 55001); Lf = interp1(lamPM, LamPM(:,8), lf, 'spline') - P.Lam;
zc = find(sign(Lf(1:end-1)) ~= sign(Lf(2:end)));
fprintf('[4] LP09 (nu = 15) is phase matched at Lambda = %d um at: %s nm\n', P.Lam, sprintf('%.0f ', lf(zc)*1e3));
write_csv(fullfile(outdir,'dual_resonance_nu15.csv'), 'lambda_nm', lf(zc(:)).'*1e3);
if DO_PLOT
figure('visible','off'); hold on;
for m = 1:8, plot(LamPM(:,m), lamPM*1e3, 'LineWidth',1.3, 'Color', C8(m,:)); end
plot([P.Lam P.Lam], [600 1800], 'k--');
plot(tap(:,1), tap(:,2), 'ko', 'MarkerSize', 7, 'LineWidth', 1.2);
xlabel('Grating period \Lambda (\mum)'); ylabel('Resonant wavelength (nm)');
xlim([100 800]); ylim([600 1800]); grid on; tt('Phase-matching curves'); box on;
legend([arrayfun(mlabel, 1:8, 'UniformOutput', false), {sprintf('\\Lambda = %d \\mum', P.Lam), 'turning point'}], 'Location','northwest');
print(fullfile(outdir,'fig4_phase_matching.png'),'-dpng','-r300');
% alternative orientation used in part of the literature (period vs wavelength)
figure('visible','off'); hold on;
for m = 1:8, plot(lamPM*1e3, LamPM(:,m), 'LineWidth',1.3, 'Color', C8(m,:)); end
plot([600 1800], [P.Lam P.Lam], 'k--');
plot(tap(:,2), tap(:,1), 'ko', 'MarkerSize', 7, 'LineWidth', 1.2);
ylabel('Grating period \Lambda (\mum)'); xlabel('Resonant wavelength (nm)');
ylim([100 800]); xlim([600 1800]); grid on; tt('Phase-matching curves'); box on;
legend([arrayfun(mlabel, 1:8, 'UniformOutput', false), {sprintf('\\Lambda = %d \\mum', P.Lam), 'turning point'}], 'Location','northwest');
print(fullfile(outdir,'fig4b_phase_matching_swapped.png'),'-dpng','-r300');
end

%% ---- 5) Transmission spectrum (Fig. 5, Table 2) -------------------------
[lam, T, res] = lpg_spectrum(P, M, 0.78, 1.22);
fprintf('[5] resonances at Lambda = %d um (nu, lambda_r [nm], kappa [1/um], kappa*L [rad], Tmin [dB], FWHM [nm], FWHM estimate [nm])\n', P.Lam);
for j = 1:size(res,1)
    fprintf('      nu=%2d  %8.1f  %.3e  %5.2f  %7.2f  %6.1f  %6.1f\n', res(j,:));
end
write_csv(fullfile(outdir,'resonances.csv'), 'nu,lambda_r_nm,kappa_per_um,kappaL_rad,Tmin_dB,FWHM_nm,FWHM_est_nm', res);
write_csv(fullfile(outdir,'spectrum.csv'), 'lambda_nm,T,T_dB', [lam(:)*1e3 T(:) 10*log10(max(T(:),1e-6))]);
ttl = sprintf('LPG transmission: a_1 = %.2f \\mum, \\Lambda = %d \\mum, L = %g mm, n_3 = %.2f', P.a1, P.Lam, P.L/1e3, P.n3);
if DO_PLOT
figure('visible','off'); plot(lam*1e3, 10*log10(max(T,1e-6)), 'LineWidth',1.2); hold on;
for j = 1:size(res,1)
    if res(j,2) >= 800 && res(j,2) <= 1200
        ylab = 2.5; if j > 1 && res(j,2) - res(j-1,2) < 12, ylab = 5.5; end   % stagger close labels
        text(max(res(j,2), 804), ylab, sprintf('%d', res(j,1)), 'HorizontalAlignment','center', 'FontSize', 9);
    end
end
text(1070, 2.5, 'numbers: mode \nu', 'HorizontalAlignment','center', 'FontSize', 9);
xlim([800 1200]); ylim([-60 8]); xlabel('Wavelength (nm)'); ylabel('Transmission (dB)');
tt(ttl); grid on; print(fullfile(outdir,'fig5_spectrum_dB.png'),'-dpng','-r300');
figure('visible','off'); plot(lam*1e3, 100*T, 'LineWidth',1.2);
xlim([800 1200]); ylim([0 105]); xlabel('Wavelength (nm)'); ylabel('Transmission (%)');
tt(ttl); grid on; print(fullfile(outdir,'fig5b_spectrum_percent.png'),'-dpng','-r300');
end

%% ---- 6) SRI sweep (Figs. 6-7, Table 5) ----------------------------------
sri   = P.sri; modes = P.sri_modes;
lamS  = linspace(0.78, 1.30, 53);
lrall = sri_sweep(P, [sri P.sri_ext], modes, lamS);   % last row = extended point
lr    = lrall(1:numel(sri), :);
for s = 1:numel(sri)
    fprintf('[6] SRI = %.2f : lambda_r (nm) = %s\n', sri(s), sprintf('%9.2f', lr(s,:)*1e3));
end
fprintf('[6] extended point SRI = %.2f : lambda_r (nm) = %s  (not used in the main analysis)\n', P.sri_ext, sprintf('%9.2f', lrall(end,:)*1e3));
shift = (lr - repmat(lr(1,:), numel(sri), 1))*1e3;              % nm, relative to air
sens  = diff(lr*1e3)./repmat(diff(sri(:)), 1, numel(modes));     % nm/RIU (finite difference)
smid  = 0.5*(sri(1:end-1)' + sri(2:end)');
% operating range: SRI where |sensitivity| >= 10 nm/RIU (log-linear interpolation)
Sthr = 10; oprange = nan(numel(modes), 3);
for j = 1:numel(modes)
    a = abs(sens(:,j)); i1 = find(a >= Sthr, 1, 'first');
    if isempty(i1), continue; end
    if i1 == 1, lo = sri(1);
    else, lo = interp1(log10(a([i1-1 i1])), smid([i1-1 i1]), log10(Sthr)); end
    oprange(j,:) = [2*modes(j)-1, lo, sri(find(~isnan(lr(:,j)),1,'last'))];
end
hdr = 'SRI,nu9,nu11,nu13,nu15';
write_csv(fullfile(outdir,'sri_sweep.csv'), hdr, [sri(:) lr*1e3]);
write_csv(fullfile(outdir,'sri_shift.csv'), hdr, [sri(:) shift]);
write_csv(fullfile(outdir,'sri_sensitivity.csv'), 'SRI_mid,nu9,nu11,nu13,nu15', [smid sens]);
write_csv(fullfile(outdir,'sri_extended.csv'), 'SRI,lr_nu9_nm,lr_nu11_nm,lr_nu13_nm,lr_nu15_nm,shift_nu9,shift_nu11,shift_nu13,shift_nu15', ...
          [P.sri_ext lrall(end,:)*1e3 (lrall(end,:)-lr(1,:))*1e3]);
write_csv(fullfile(outdir,'operating_range.csv'), 'nu,SRI_low_10nm_per_RIU,SRI_high', oprange);
fprintf('[6] shift at SRI %.2f (nm): %s\n', sri(end), sprintf('%8.2f', shift(end,:)));
fprintf('[6] |S| >= %g nm/RIU: %s\n', Sthr, sprintf('nu=%d: %.3f-%.2f; ', oprange.'));
fprintf('[6] |S| ratio (last interval %.2f-%.2f / first interval): %s\n', sri(end-1), sri(end), sprintf('%.0f ', abs(sens(end,:)./sens(1,:))));
% analytical scaling of the sensitivity with the surrounding index:
%   S ~ n3*(n2^2 - n3^2)^(-3/2)  (derivative of the shift formula of Chiang et al., 2000)
g  = @(x) x.*(P.n2^2 - x.^2).^(-1.5);  g0 = @(x) (P.n2^2 - x.^2).^(-1.5);
fprintf('[6] predicted ratio: n3*(n2^2-n3^2)^(-3/2) -> %.0f ; (n2^2-n3^2)^(-3/2) alone -> %.0f\n', ...
        g(smid(end))/g(smid(1)), g0(smid(end))/g0(smid(1)));
if DO_PLOT
srip = [1.00 1.33 1.40 1.44]; Tq = cell(1, numel(srip));
for s = 1:numel(srip), Q = P; Q.n3 = srip(s); [lq, Tq{s}] = lpg_spectrum(Q, M, 0.78, 1.22); end
lg = arrayfun(@(v) sprintf('SRI = %.2f', v), srip, 'UniformOutput', false);
figure('visible','off'); hold on;
for s = 1:numel(srip), plot(lq*1e3, 10*log10(max(Tq{s},1e-6)), 'LineWidth', 1.1, 'Color', C8(s,:)); end
xlim([800 1200]); ylim([-60 5]); xlabel('Wavelength (nm)'); ylabel('Transmission (dB)'); grid on; box on;
legend(lg, 'Location','southwest'); tt('Transmission spectra for different SRI');
print(fullfile(outdir,'fig6_sri_spectra_dB.png'),'-dpng','-r300');
figure('visible','off'); hold on;
for s = 1:numel(srip), plot(lq*1e3, 100*Tq{s}, 'LineWidth', 1.1, 'Color', C8(s,:)); end
xlim([800 1200]); ylim([0 105]); xlabel('Wavelength (nm)'); ylabel('Transmission (%)'); grid on; box on;
legend(lg, 'Location','southwest'); tt('Transmission spectra for different SRI');
print(fullfile(outdir,'fig6b_sri_spectra_percent.png'),'-dpng','-r300');
mk = {'-o','-s','-^','-d'};
figure('visible','off'); hold on;
for j = 1:numel(modes), plot(sri, shift(:,j), mk{j}, 'LineWidth',1.3, 'Color', C8(modes(j),:)); end
xlim([1.0 1.45]); xlabel('Surrounding refractive index'); ylabel('Resonance wavelength shift (nm)'); grid on; box on;
legend(arrayfun(mlabel, modes, 'UniformOutput', false), 'Location','southwest');
tt('Resonance wavelength shift vs. SRI'); print(fullfile(outdir,'fig7_sri_shift.png'),'-dpng','-r300');
figure('visible','off');
for j = 1:numel(modes), semilogy(smid, abs(sens(:,j)), mk{j}, 'LineWidth',1.3, 'Color', C8(modes(j),:)); hold on; end
semilogy([1 1.45], [Sthr Sthr], 'k--'); xlim([1.0 1.45]);
xlabel('Surrounding refractive index'); ylabel('|d\lambda_r/dn_3| (nm/RIU)'); grid on;
legend([arrayfun(mlabel, modes, 'UniformOutput', false), {'10 nm/RIU'}], 'Location','northwest');
tt('SRI sensitivity'); print(fullfile(outdir,'fig7b_sri_sensitivity.png'),'-dpng','-r300');
end
fprintf('lpg_main done in %.1f s. Figures and CSV files are in ./%s\n', toc(t_start), outdir);
