% FOM_MAP  Sensitivity, bandwidth and figure of merit versus grating period.
%
% Purpose (paper, Results Section 2.3): show that proximity to
% the dispersion turning point (TAP) raises the SRI sensitivity S and the
% bandwidth FWHM by the same factor 1/|dn_g|, so that their ratio
%       FOM = |S|/FWHM  ~=  pi*L*|dncl/dn3| / (2*x_h(kappa*L)*lambda_r)      Eq. (7)
% does not depend on how close the grating is to the TAP. dn_g is the
% group-index difference of the core and cladding modes, dncl/dn3 the
% derivative of the cladding-mode effective index with respect to the SRI,
% and x_h the half-width (in delta*L) of a band with coupling strength kappa*L.
%
% For every cladding mode the phase-matching curve is parametrised by the
% resonant wavelength: each lambda on the grid corresponds to the period
% Lambda = lambda/(nco - ncl) that puts the resonance there. For each such
% design the script
%   - computes S = d lambda_r / d n3 exactly by re-solving phase matching at
%     n3 +/- h (central difference),
%   - simulates the band with eq. (4) on a fine grid and measures FWHM at
%     half of the coupled power (same definition as band_fwhm.m),
%   - evaluates the closed form, Eq. (7), and the first-order estimates
%     S1 = -lambda*dncl/dn3/dn_g and FWHM1 = 2*x_h*lambda^2/(pi*L*|dn_g|).
% Output: results/fom_map.csv (all designs). Figure 8 of the paper: python make_fig8.py
% Run time: a few minutes in MATLAB.
clear; P = lpg_params(); outdir = 'results'; if ~exist(outdir,'dir'), mkdir(outdir); end
n3  = 1.33;                 % water (design point); other n3 via n3list below
h   = 0.005;                % step of the central difference in n3
modes = 3:9;                % cladding modes j: LP04..LP0,10 = nu 5,7,...,17
lc  = linspace(0.70, 1.80, 111);   % coarse grid of resonant wavelengths (um)
nr  = 6001;
n3list = [1.00 1.33 1.40];
out = zeros(0, 14);
for n3 = n3list
  NCO = zeros(numel(lc),1); NCL = nan(numel(lc), max(modes), 3); KAP = nan(numel(lc), max(modes));
  for i = 1:numel(lc)
    [n1, n2] = index_model(P, lc(i));
    NCO(i) = neff_core(lc(i), P.a1, n1, n2);
    for s = 1:3
      NCL(i,:,s) = neff_clad(lc(i), P.a1, P.a2, n1, n2, n3 + (s-2)*h, max(modes));
    end
    KAP(i,:) = coupling(lc(i), NCO(i), NCL(i,:,2), P.a1, P.a2, n1, n2, n3, P.dn, nr);
  end
  for j = modes
    ok = all(~isnan(squeeze(NCL(:,j,:))), 2) & ~isnan(KAP(:,j));
    l = lc(ok); if numel(l) < 8, continue; end
    dn0 = NCO(ok) - NCL(ok,j,2); dnm = NCO(ok) - NCL(ok,j,1); dnp = NCO(ok) - NCL(ok,j,3);
    kap = KAP(ok,j);
    % dense spline representation of the phase-matching functions
    lf  = linspace(l(1), l(end), 4001).';
    D0  = interp1(l, dn0, lf, 'spline'); Dm = interp1(l, dnm, lf, 'spline'); Dp = interp1(l, dnp, lf, 'spline');
    K   = interp1(l, kap, lf, 'spline');
    dng = D0 - lf.*gradient(D0, lf);
    dcl = -(Dp - Dm)/(2*h);                  % d ncl / d n3 (= -d dn / d n3)
    Lam = lf./D0;                            % period that puts the resonance at lf
    for q = 20:40:numel(lf)-20               % ~100 designs per mode
      LamQ = Lam(q); lr = lf(q);
      % exact sensitivity: resonance nearest lr at n3 -/+ h with the same period
      fm = Dm.*LamQ - lf; fp = Dp.*LamQ - lf;
      rm = zcross(lf, fm, lr); rp = zcross(lf, fp, lr);
      S = (rp - rm)/(2*h)*1e3;               % nm/RIU
      % simulated band (eq. 4) and its FWHM at half of the coupled power
      delta = pi*(D0./lf - 1/LamQ);
      ss = sqrt(K.^2 + delta.^2);
      Tm = 1 - (K.^2./ss.^2).*sin(ss*P.L).^2;
      [~, qq] = min(abs(lf - lr));
      FW = band_fwhm(lf.', Tm.', qq);       % nm
      kL = K(q)*P.L; xh = xhalf(kL);
      S1  = -lr*dcl(q)/dng(q)*1e3;
      FW1 = 2*xh*lr^2/(pi*P.L*abs(dng(q)))*1e3;
      FOMA = pi*P.L*abs(dcl(q))/(2*xh*lr);   % closed form, Eq. (7), 1/RIU
      out(end+1,:) = [n3, 2*j-1, LamQ, lr*1e3, dng(q), dcl(q), kL, S, FW, S1, FW1, abs(S)/FW, FOMA, 10*log10(max(Tm(qq),1e-6))]; %#ok<AGROW>
    end
  end
  fprintf('n3 = %.2f done\n', n3);
end
write_csv(fullfile(outdir,'fom_map.csv'), ...
  'n3,nu,Lambda_um,lambda_r_nm,dn_g,dncl_dn3,kappaL,S_nm_RIU,FWHM_nm,S_firstorder,FWHM_firstorder,FOM_sim,FOM_closedform,Tmin_dB', out);

% ---- figure 8 (MATLAB). Without a graphics toolkit (headless GNU Octave) use
%      python make_fig8.py, which draws Figure 8 of the paper from results/fom_map.csv
try
sel = out(:,1) == 1.33 & ~isnan(out(:,9));
nus = unique(out(sel,2)).'; cols = lines(numel(nus));
figure('visible','off','position',[0 0 1500 450]);
lab = {'|S| (nm/RIU)', 'FWHM (nm)', 'FOM = |S|/FWHM (RIU^{-1})'}; colx = [8 9 12];
for p = 1:3
  subplot(1,3,p); hold on;
  for k = 1:numel(nus)
    r = sel & out(:,2) == nus(k);
    semilogy(abs(out(r,5)), abs(out(r,colx(p))), 'o-', 'color', cols(k,:), 'markersize', 3);
  end
  set(gca,'yscale','log','xscale','log'); grid on; box on;
  xlabel('|\Deltan_g| (group-index difference)'); ylabel(lab{p});
  if p == 1, legend(arrayfun(@(v) sprintf('\\nu = %d', v), nus, 'uniformoutput', false), 'location', 'northeast'); end
end
print(fullfile(outdir,'fig8_quicklook_matlab.png'), '-dpng', '-r200');   % paper figure: python make_fig8.py
catch err
  fprintf('fom_map: figure skipped (%s); run python make_fig8.py\n', err.message);
end
fprintf('fom_map: %d designs -> results/fom_map.csv (figure 8: python make_fig8.py)\n', size(out,1));

