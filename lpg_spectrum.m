function [lam, T, res, Tm_all] = lpg_spectrum(P, M, lmin, lmax)
% Transmission spectrum, eqs. (2) and (4).
% Each cladding mode is treated as an independent two-mode (core-cladding)
% coupling problem and the total transmission is the product over the M
% modes. This is valid when the attenuation bands do not overlap.
% Mode indices and coupling coefficients are computed on a coarse wavelength
% grid (P.nc points, default 45) and interpolated with cubic splines onto a
% fine grid (P.nf points, default 40,001). The overlap integrals use P.nr
% radial points (default 12,001). convergence_check.m doubles all three.
%
% res columns: nu, lambda_r (nm, zero of the detuning by linear interpolation
%   between fine-grid points), kappa (1/um), kappa*L (rad), Tmin (dB, value
%   of eq. (4) at the grid point nearest to phase matching),
%   FWHM (nm): width at half of the coupled power (= 3-dB bandwidth for a
%              fully coupled band); NaN if the band is cut by the window
%   FWHM_est (nm): analytical estimate 0.8*lambda^2/(L*|dn_g|) commonly used
%              for LPGs (dn_g = group-index difference; usually attributed
%              to Vengsarkar et al., J. Lightw. Technol. 14, 58 (1996) -
%              check the citation before using it in the paper). It is a
%              cross-check only and diverges at the turning point
dndc = pget(P, 'dn_dc', 0);
nc_grid = pget(P, 'nc', 45); nf = pget(P, 'nf', 40001); nr = pget(P, 'nr', 12001);
lc  = linspace(lmin, lmax, nc_grid);
nco = zeros(size(lc)); ncl = zeros(numel(lc), M); kap = ncl;
for i = 1:numel(lc)
    [n1, n2] = index_model(P, lc(i));
    nco(i)   = neff_core(lc(i), P.a1, n1, n2);
    ncl(i,:) = neff_clad(lc(i), P.a1, P.a2, n1, n2, P.n3, M);
    kap(i,:) = coupling(lc(i), nco(i), ncl(i,:), P.a1, P.a2, n1, n2, P.n3, P.dn, nr);
    if dndc ~= 0      % DC self-coupling (optional, P.dn_dc; default 0 = model of the paper)
        [Gco, Gcl] = confinement(lc(i), nco(i), ncl(i,:), P.a1, P.a2, n1, n2, P.n3, nr);
        nco(i) = nco(i) + dndc*Gco; ncl(i,:) = ncl(i,:) + dndc*Gcl;
    end
end
lam = linspace(lmin, lmax, nf);
T = ones(size(lam)); res = zeros(0, 7); Tm_all = ones(M, nf);
nc = interp1(lc, nco, lam, 'spline');
for m = 1:M
    if any(isnan(ncl(:,m))), continue; end
    nl = interp1(lc, ncl(:,m), lam, 'spline');
    k  = interp1(lc, kap(:,m), lam, 'spline');
    delta = pi*(nc - nl)./lam - pi/P.Lam;              % detuning, eq. (2)
    s = sqrt(k.^2 + delta.^2);
    Tm = 1 - (k.^2./s.^2).*sin(s*P.L).^2;               % eq. (4)
    T = T.*Tm; Tm_all(m,:) = Tm;
    dneff = nc - nl;
    dng = dneff - lam.*gradient(dneff, lam);            % group-index difference
    z = find(sign(delta(1:end-1)) ~= sign(delta(2:end)));
    for q = z(:).'
        lr = lam(q) - delta(q)*(lam(q+1)-lam(q))/(delta(q+1)-delta(q));   % delta = 0
        res(end+1,:) = [2*m-1, lr*1e3, k(q), k(q)*P.L, ...
                        10*log10(max(Tm(q),1e-6)), band_fwhm(lam, Tm, q), ...
                        0.8*lam(q)^2/(P.L*abs(dng(q)))*1e3]; %#ok<AGROW>
    end
end
if ~isempty(res), res = sortrows(res, 2); end
end
