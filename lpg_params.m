function P = lpg_params()
% LPG_PARAMS  Single source of the model parameters (Table 1 of the paper).
% All scripts call this function so that every table and figure uses the
% same parameter set. Units: micrometres.
%
% How the fibre parameters were chosen (state this in the paper):
%   a1, n1 and n2 were CHOSEN so that NA = sqrt(n1^2 - n2^2) = 0.121 and the
%   cut-off wavelength (646 nm) fall inside the PS750 specification
%   (NA 0.12-0.14, cut-off 610-750 nm). They are not values printed in the
%   datasheet. The datasheet comparison in Table 3 is therefore a
%   consistency check, not an independent validation.
P.a1   = 2.05;     % core radius
P.a2   = 62.5;     % cladding radius
P.n1   = 1.4583;   % core index (at the reference wavelength when disp = false)
P.n2   = 1.4533;   % cladding index (= fused silica at ~0.80 um)
P.n3   = 1.0000;   % surrounding index (air)
P.dn   = 2e-4;     % AMPLITUDE of the index modulation, uniform over the core
P.L    = 30e3;     % grating length (30 mm)
P.Lam  = 400;      % grating period
P.M    = 20;       % number of cladding modes: LP02..LP0,21 (Erdogan nu = 1,3,...,39)
P.lam0 = 1.0;      % reference wavelength for Figs. 1-3
P.disp = false;    % true: silica Sellmeier dispersion, constant n1-n2 step
% Surrounding refractive index (SRI) values of the main analysis: air to 1.44,
% the range of Patrick et al., J. Lightw. Technol. 16, 1606 (1998). SRI = 1.45
% is kept only as an extended point: with material dispersion it exceeds the
% cladding index at the nu = 15 band (leaky regime), see dispersion_check.m.
P.sri       = [1.00:0.05:1.40, 1.41:0.01:1.44];   % 13 values
P.sri_ext   = 1.45;                               % reported separately
P.sri_modes = [5 6 7 8];                          % 5th-8th cladding modes: LP06..LP09 = nu 9, 11, 13, 15
% Manufacturer specification of PS750 used in Table 3 (verify with datasheet [19])
P.spec_MFD780 = [4.4 5.9];   % um, mode-field diameter at 780 nm
P.spec_cutoff = [610 750];   % nm
P.spec_NA     = [0.12 0.14];
end
